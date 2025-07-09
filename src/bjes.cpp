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

    const std::string output_ntuple_name = output_folder_ + "/ntuples/ntuples_"
                                                        + std::to_string(sample_label) + ".root";

    std::vector<std::string> output_variables = {"ljet_pt", "ljet_eta", "ljet_phi", "ljet_m", "ljet_e", "n_ljets",
                                                 "jet_pt", "jet_eta", "jet_phi", "jet_m", "jet_e", "jet_jvt",
                                                 "ljet_bJR10v00_mass", "ljet_bJR10v00_pt",
                                                 "ljet_bJR10v01_mass", "ljet_bJR10v01_pt",
                                                 "dhbb", "n_jets"};

    for (auto wp : config_.analysis.wps) {
      std::string pass_dhbb = "pass_" + wp;
      output_variables.emplace_back(pass_dhbb);
    }


    if (!is_data) {
      output_variables.insert(output_variables.end(), {"ljet_truth_label", "jet_parton_truth_label", "jet_eff_jvt"});
    }

    if (config_.analysis.analysis == "zbby") {
      output_variables.insert(output_variables.end(), {"ph_pt", "ph_eta", "ph_phi", "met_met", "met_phi"});
    }
    df_.value().Snapshot("tree", output_ntuple_name, output_variables);
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

  void Analysis::define_physics_variables(bool is_data)
  {
    const std::string ljet_prefix = "recojet_antikt10UFO_";
    const std::string jet_prefix = "recojet_antikt4PFlow_";
    const std::string gn2x_prefix = ljet_prefix + "GN2Xv01_";

    std::vector<std::pair<std::string, std::string>> variables = {
      {"ljet_m", ljet_prefix + "m_NOSYS/1000"},
      {"ljet_pt", ljet_prefix + "pt_NOSYS/1000"},
      {"ljet_phi", ljet_prefix + "phi_NOSYS"},
      {"ljet_eta", ljet_prefix + "eta_NOSYS"},

      {"jet_m", jet_prefix + "m_NOSYS/1000"},
      {"jet_pt", jet_prefix + "pt_NOSYS/1000"},
      {"jet_phi", jet_prefix + "phi_NOSYS"},
      {"jet_eta", jet_prefix + "eta_NOSYS"},
      {"jet_jvt", jet_prefix + "Jvt_NOSYS"},

      {"ljet_bJR10v00_mass", ljet_prefix + "bJR10v00_mass_NOSYS/1000"},
      {"ljet_bJR10v01_mass", ljet_prefix + "bJR10v01_mass_NOSYS/1000"},
      {"ljet_bJR10v00_pt", ljet_prefix + "bJR10v00_pt_NOSYS/1000"},
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

    std::vector<std::pair<std::string, std::string>> mc_variables = {
      {"ljet_truth_label", ljet_prefix + "R10TruthLabel_R22v1_NOSYS"},
      {"jet_parton_truth_label", jet_prefix + "PartonTruthLabelID_NOSYS"},
      {"jet_eff_jvt", jet_prefix + "jvt_effSF_NOSYS"},
    };

    for (const auto& [name, expr] : variables) {
      df_ = df_.value().Define(name, expr);
    }

    if (!is_data) {
      for (const auto& [name, expr] : mc_variables) {
        df_ = df_.value().Define(name, expr);
      }
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
  }
}
