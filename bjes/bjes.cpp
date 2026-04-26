#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <map>
#include <utility>
#include <algorithm>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TF1.h>
#include <TChain.h>
#include <spdlog/spdlog.h>
#include "output_manager.h"
#include "input_handler.h"
#include "weights.h"
#include "config_bjes.h"
#include "bjes.h"
#include "gn2x.h"

namespace bjes {
  Analysis::Analysis(const config_bjes::Config& config,
                     const std::string& input_folder,
                     std::string output_folder,
                     int n_threads)
    : config_(config),
      input_folder_(input_folder),
      output_folder_(output_folder),
      n_threads_(n_threads)
  {
  }

  void Analysis::run()
  {
  output_folder_ = output_manager::create_output_folder(output_folder_);
    output_manager::set_logger(output_folder_);

    spdlog::info("Running analysis...");
    if (n_threads_ > 1) {
      ROOT::EnableImplicitMT(n_threads_);
      spdlog::info("Multi-threading enabled using: {} threads.", n_threads_);
    }
    else {
      spdlog::info("Multi-threading disabled, RDataFrame will run on a single core.");
    }
    spdlog::info("Running on MC samples...");
    for (auto dsid : config_.analysis.dsids) {
      run_sample(dsid, false);
    }
    spdlog::info("Running on data samples...");
    for (auto year : config_.analysis.years) {
      run_sample(year, true);
    }
    spdlog::info("Analysis completed successfully.");
  }

  void Analysis::run_sample(int sample_label, bool is_data)
  {
    if (is_data) spdlog::info("Creating RDataFrame for data sample for year: {}", sample_label);
    else         spdlog::info("Creating RDataFrame for MC sample with dsid: {}",  sample_label);

    auto input_handler = InputHandler(input_folder_, config_.ntuple.tree_name, sample_label, is_data, config_.analysis.metadata);
    auto chain_ptr = input_handler.get_chain();

    ROOT::RDataFrame df(*chain_ptr);
    ROOT::RDF::Experimental::AddProgressBar(df);
    df_ = df;

    update_cutflow("Initial");
    apply_triggers();
    define_physics_variables(is_data);
    GN2XHandler gn2x_handler_(config_.analysis.flatmass);
    df_ = gn2x_handler_.define_pass_rvec(df_.value(), "dhbb", "ljet_m", "ljet_pt");

    const std::string output_file_name = output_folder_ + "/histograms/hists_"
                                                        + std::to_string(sample_label);

    const std::string output_file_name_full = output_file_name + ".root";
    const std::string output_file_name_mc20a = output_file_name + "_mc20a.root";
    const std::string output_file_name_mc20d = output_file_name + "_mc20d.root";
    const std::string output_file_name_mc20e = output_file_name + "_mc20e.root";
    // mc23
    const std::string output_file_name_mc23a = output_file_name + "_mc23a.root";
    const std::string output_file_name_mc23d = output_file_name + "_mc23d.root";
    const std::string output_file_name_mc23e = output_file_name + "_mc23e.root";

    const std::string output_ntuple_name = output_folder_ + "/ntuples/ntuples_"
                                                        + std::to_string(sample_label) + ".root";

    std::vector<std::string> output_variables = {"ljet_pt", "ljet_eta", "ljet_phi", "ljet_m", "ljet_e", "n_ljets",
                                                 "jet_pt", "jet_eta", "jet_phi", "jet_m", "jet_e", "jet_jvt",
                                                 "ljet_bJR10v00_mass", "ljet_bJR10v00_pt",
                                                 "ljet_bJR10v01_mass", "ljet_bJR10v01_pt",
                                                 "dhbb", "n_jets", "total_weight"};

    for (auto wp : config_.analysis.wps) {
      std::string pass_dhbb = "pass_" + wp;
      output_variables.emplace_back(pass_dhbb);
    }


    if (!is_data) {
      output_variables.insert(output_variables.end(), {
      "ljet_truth_label", "jet_parton_truth_label", "jet_eff_jvt",
      // --- NEW ---
      "matched_tjet_pt", "matched_tjet_m",
      "matched_tjet_eta", "matched_tjet_phi",
      "matched_tjet_dR",
      "ljet_n_bhadrons",
      "ljet_n_chadrons"
    });
  }

    if (config_.analysis.analysis == "zbby") {
      output_variables.insert(output_variables.end(), {"ph_pt", "ph_eta", "ph_phi", "ph_e", "met_met", "met_phi"});
    }
    df_.value().Snapshot("nominal", output_ntuple_name, output_variables);
  }

  void Analysis::update_cutflow(std::string label)
  {
    auto sum_weights = df_.value().Sum("total_weight");
    auto n_entries = df_.value().Count();
    std::pair<std::string, double> cut (label, sum_weights.GetValue());
    std::pair<std::string, double> cut_unweighted (label, n_entries.GetValue());
    cutflow_.emplace_back(cut);
    cutflow_unweighted_.emplace_back(cut_unweighted);
  }

  void Analysis::apply_triggers()
  {
    // Build trigger condition expression
    std::string trigger_expr;

    for (const auto& [year, triggers] : config_.analysis.trigger_map) {
      std::string condition = "(";
      for (size_t i = 0; i < triggers.size(); ++i) {
        std::string trigger_name = "trigPassed_" + triggers[i];
        df_ = df_.value().DefaultValueFor(trigger_name, false);
        condition += trigger_name;
        if (i < triggers.size() - 1) condition += " || ";
      }
      condition += ")";
      trigger_expr += "(dataTakingYear == " + std::to_string(year) + " && " + condition + ") || ";
    }

    // Remove the extra ||
    trigger_expr = trigger_expr.substr(0, trigger_expr.size() - 4);

    spdlog::info("Appplying trigger selection : {}", trigger_expr);
    df_ = df_.value().Filter(trigger_expr);

    update_cutflow("Triggers");
  }

  std::string Analysis::get_branch_name(const std::string& base_name,
                                        const std::vector<std::string>& candidates)
  {
    for (const auto& candidate : candidates) {
      if (df_.value().HasColumn(candidate)) {
        spdlog::info("Using branch '{}' for '{}'", candidate, base_name);
        return candidate;
      }
    }
    // Build the list of tried branches manually
    std::string tried;
    for (size_t i = 0; i < candidates.size(); ++i) {
      tried += candidates[i];
      if (i < candidates.size() - 1) tried += ", ";
    }
    spdlog::error("No valid branch found for '{}'. Tried: {}", base_name, tried);
    throw std::runtime_error("Missing required branch: " + base_name);
  }

  void Analysis::define_physics_variables(bool is_data)
  {
    const std::string ljet_prefix = "recojet_antikt10UFO_";
    const std::string jet_prefix = "recojet_antikt4PFlow_";
    const std::string gn2x_prefix = ljet_prefix + "GN2Xv01_";

    // Resolve branch name variants
    auto bjr10v00_mass = get_branch_name("ljet_bJR10v00_mass", {
      ljet_prefix + "bJR10v00Ext_mass_NOSYS",
      ljet_prefix + "bJR10v00_mass_NOSYS"
    });
    auto bjr10v00_pt = get_branch_name("ljet_bJR10v00_pt", {
      ljet_prefix + "bJR10v00Ext_pt_NOSYS",
      ljet_prefix + "bJR10v00_pt_NOSYS"
    });

    // JVT branch resolution:
    auto jet_jvt_branch = get_branch_name("jet_jvt", {
        jet_prefix + "jvt_selection_NOSYS",   // zbby
        jet_prefix + "Jvt_NOSYS"              // zbbj
    });

    std::vector<std::pair<std::string, std::string>> variables = {
      {"ljet_m", ljet_prefix + "m_NOSYS/1000"},
      {"ljet_pt", ljet_prefix + "pt_NOSYS/1000"},
      {"ljet_phi", ljet_prefix + "phi_NOSYS"},
      {"ljet_eta", ljet_prefix + "eta_NOSYS"},

      {"jet_m", jet_prefix + "m_NOSYS/1000"},
      {"jet_pt", jet_prefix + "pt_NOSYS/1000"},
      {"jet_phi", jet_prefix + "phi_NOSYS"},
      {"jet_eta", jet_prefix + "eta_NOSYS"},
      {"jet_jvt", jet_jvt_branch},
      //{"recojet_antikt4PFlow_jvt_NOSYS", jet_jvt_branch},
      //recojet_antikt4PFlow_jvt_selection_NOSYS in Zbby and recojet_antikt4PFlow_jvt_NOSYS in Zbbj
      //{"jet_jvt", jet_prefix + "Jvt_NOSYS"}, 

      //{"ljet_bJR10v00_mass", ljet_prefix + "bJR10v00Ext_mass_NOSYS/1000"}, //bJR10v00Ext_mass_NOSYS or bJR10v00_mass_NOSYS
      //{"ljet_bJR10v01_mass", ljet_prefix + "bJR10v01_mass_NOSYS/1000"},
      //{"ljet_bJR10v00_pt", ljet_prefix + "bJR10v00Ext_pt_NOSYS/1000"},
      //{"ljet_bJR10v01_pt", ljet_prefix + "bJR10v01_pt_NOSYS/1000"},

      {"ljet_bJR10v00_mass", bjr10v00_mass + "/1000"},
      {"ljet_bJR10v00_pt", bjr10v00_pt + "/1000"},
      {"ljet_bJR10v01_mass", ljet_prefix + "bJR10v01_mass_NOSYS/1000"},
      {"ljet_bJR10v01_pt", ljet_prefix + "bJR10v01_pt_NOSYS/1000"},

      {"ljet_phbb", gn2x_prefix + "phbb_NOSYS"},
      {"ljet_phcc", gn2x_prefix + "phcc_NOSYS"},
      {"ljet_ptop", gn2x_prefix + "ptop_NOSYS"},
      {"ljet_pqcd", gn2x_prefix + "pqcd_NOSYS"},

      {"n_ljets", "ljet_pt.size()"},
      {"n_jets", "jet_pt.size()"},
    };

    std::vector<std::pair<std::string, std::string>> zbby_variables = {
      {"ph_pt", "ph_pt_NOSYS/1000"},
      {"ph_eta", "ph_eta_NOSYS"},
      {"ph_phi", "ph_phi_NOSYS"},

      {"met_met", "met_met_NOSYS"},
      {"met_phi", "met_phi_NOSYS"},
    };

    for (const auto& [name, expr] : variables) {
      df_ = df_.value().Define(name, expr);
    }

    if (!is_data) {
      // JVT EFF branch resolution, MC only:
      auto jet_eff_jvt_branch = get_branch_name("jet_eff_jvt", {
          jet_prefix + "jvt_effSF_NOSYS",
          "jvt_effSF_NOSYS"
      });
      
      std::vector<std::pair<std::string, std::string>> mc_variables = {
        {"ljet_truth_label", ljet_prefix + "R10TruthLabel_R22v1_NOSYS"},
        {"jet_parton_truth_label", jet_prefix + "PartonTruthLabelID_NOSYS"},
        //{"jet_eff_jvt", jet_prefix + "jvt_effSF_NOSYS"},
        {"jet_eff_jvt", jet_eff_jvt_branch},
        {"ljet_n_bhadrons",     ljet_prefix + "GhostBHadronsFinalCount_NOSYS"},  
        {"ljet_n_chadrons",     ljet_prefix + "GhostCHadronsFinalCount_NOSYS"},  
      };

      for (const auto& [name, expr] : mc_variables) {
        df_ = df_.value().Define(name, expr);
      }
      
        const std::string tjet_prefix = "truthjet_antikt10SoftDrop_";

    // alias raw truth branches to short names
    /*df_ = df_.value()
      .Alias("tjet_eta", "truthjet_antikt10SoftDrop_eta")
      .Alias("tjet_phi", "truthjet_antikt10SoftDrop_phi")
      .Define("tjet_pt", [](const ROOT::RVecF& v){ return v / 1000.f; }, {"truthjet_antikt10SoftDrop_pt"})
      .Define("tjet_m",  [](const ROOT::RVecF& v){ return v / 1000.f; }, {"truthjet_antikt10SoftDrop_m"});
    */
    df_ = df_.value()
      .Define("tjet_eta", "truthjet_antikt10SoftDrop_eta")
      .Define("tjet_phi", "truthjet_antikt10SoftDrop_phi")
      .Define("tjet_pt",  "truthjet_antikt10SoftDrop_pt/1000")
      .Define("tjet_m",   "truthjet_antikt10SoftDrop_m/1000");

    // ΔR matching: for each reco ljet, find the closest truth jet
    auto get_matched = [](
        const ROOT::RVecF& reco_eta,  const ROOT::RVecF& reco_phi,
        const ROOT::RVecF& truth_eta, const ROOT::RVecF& truth_phi,
        const ROOT::RVecF& truth_vals) -> ROOT::RVecF
    {
      ROOT::RVecF result(reco_eta.size(), -999.f);
      for (size_t i = 0; i < reco_eta.size(); ++i) {
        float min_dR = std::numeric_limits<float>::max();
        int   best   = -1;
        for (size_t j = 0; j < truth_eta.size(); ++j) {
          float deta = reco_eta[i] - truth_eta[j];
          float dphi = reco_phi[i] - truth_phi[j];
          // wrap dphi into [-π, π]
          while (dphi >  M_PI) dphi -= 2.f * M_PI;
          while (dphi < -M_PI) dphi += 2.f * M_PI;
          float dR = std::sqrt(deta*deta + dphi*dphi);
          if (dR < min_dR) { min_dR = dR; best = j; }
        }
        if (best >= 0) result[i] = truth_vals[best];
      }
      return result;
    };

    // also store the minimum ΔR itself (useful for a matching quality cut later)
    auto get_min_dR = [](
        const ROOT::RVecF& reco_eta,  const ROOT::RVecF& reco_phi,
        const ROOT::RVecF& truth_eta, const ROOT::RVecF& truth_phi) -> ROOT::RVecF
    {
      ROOT::RVecF result(reco_eta.size(), -999.f);
      for (size_t i = 0; i < reco_eta.size(); ++i) {
        float min_dR = std::numeric_limits<float>::max();
        for (size_t j = 0; j < truth_eta.size(); ++j) {
          float deta = reco_eta[i] - truth_eta[j];
          float dphi = reco_phi[i] - truth_phi[j];
          while (dphi >  M_PI) dphi -= 2.f * M_PI;
          while (dphi < -M_PI) dphi += 2.f * M_PI;
          float dR = std::sqrt(deta*deta + dphi*dphi);
          if (dR < min_dR) min_dR = dR;
        }
        result[i] = (min_dR < std::numeric_limits<float>::max()) ? min_dR : -999.f;
      }
      return result;
    };

    df_ = df_.value()
      .Define("matched_tjet_pt",  get_matched, {"ljet_eta", "ljet_phi", "tjet_eta", "tjet_phi", "tjet_pt"})
      .Define("matched_tjet_m",   get_matched, {"ljet_eta", "ljet_phi", "tjet_eta", "tjet_phi", "tjet_m"})
      .Define("matched_tjet_eta", get_matched, {"ljet_eta", "ljet_phi", "tjet_eta", "tjet_phi", "tjet_eta"})
      .Define("matched_tjet_phi", get_matched, {"ljet_eta", "ljet_phi", "tjet_eta", "tjet_phi", "tjet_phi"})
      .Define("matched_tjet_dR",  get_min_dR,  {"ljet_eta", "ljet_phi", "tjet_eta", "tjet_phi"});

    }

    if (config_.analysis.analysis == "zbby") {
      for (const auto& [name, expr] : zbby_variables) {
        df_ = df_.value().Define(name, expr);
      }
    }

    auto get_energy = [](const ROOT::RVecF& pt, const ROOT::RVecF& eta, const ROOT::RVecF& m) -> ROOT::RVecF
    {
      auto cosh_eta = ROOT::VecOps::cosh(eta);
      auto p = pt * cosh_eta;
      return ROOT::VecOps::sqrt(p * p + m * m);
    };

    auto get_dhbb = [](const ROOT::RVecF& phbb, const ROOT::RVecF& phcc, const ROOT::RVecF& ptop, const ROOT::RVecF& pqcd) -> ROOT::RVecF
    {
      float fcc = 0.2;
      float ftop = 0.25;
      return ROOT::VecOps::log(phbb/(fcc*phcc + ftop*ptop + (1 - ftop - fcc) * pqcd));
    };

    df_ = df_.value()
             .Define("ljet_e", get_energy, {"ljet_pt", "ljet_eta", "ljet_m"})
             .Define("dhbb", get_dhbb, {"ljet_phbb", "ljet_phcc", "ljet_ptop", "ljet_pqcd"})
             .Define("jet_e", get_energy, {"jet_pt", "jet_eta", "jet_m"});

    if (config_.analysis.analysis == "zbby") {
      df_ = df_.value()
               .Define("ph_m", [](const ROOT::RVecF& pt){ return ROOT::RVecF(pt.size(), 0.f); }, {"ph_pt"})
               .Define("ph_e", get_energy, {"ph_pt", "ph_eta", "ph_m"});
    }
  }
}
