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
#include "config.h"
#include "xbby_calib.h"
#include "gn2x.h"

namespace xbbycalib {
  Analysis::Analysis(const config::Config& config,
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
    define_physics_variables();
    select_Z_candidate();
    get_dhbb();
    GN2XHandler gn2x_handler_(config_.analysis.flatmass);
    df_ = gn2x_handler_.define_pass(df_.value(), "zcand_gn2x_ftop0", "zcand_m", "zcand_pt");

    const std::string output_file_name = output_folder_ + "/histograms/hists_"
                                                        + std::to_string(sample_label);

    const std::string output_file_name_full = output_file_name + ".root";
    const std::string output_file_name_mc20a = output_file_name + "_mc20a.root";
    const std::string output_file_name_mc20d = output_file_name + "_mc20d.root";
    const std::string output_file_name_mc20e = output_file_name + "_mc20e.root";

    save_histograms(output_file_name_full, "total_weight");
    if (!is_data) {
      save_histograms(output_file_name_mc20a, "total_weight", "dataTakingYear == 2015 || dataTakingYear == 2016");
      save_histograms(output_file_name_mc20d, "total_weight", "dataTakingYear == 2017");
      save_histograms(output_file_name_mc20e, "total_weight", "dataTakingYear == 2018");
    }

    const std::string output_ntuple_name = output_folder_ + "/ntuples/ntuples_"
                                                        + std::to_string(sample_label) + ".root";

    std::vector<std::string> output_variables = {"zcand_m", "zcand_pt", "ph_pt", "ph_eta", "ph_phi",
                                                 "total_weight", "eventNumber",
                                                 "zcand_log_phbb", "zcand_log_phcc", "zcand_log_ptop", "zcand_log_pqcd",
                                                 "zcand_phbb", "zcand_phcc", "zcand_ptop", "zcand_pqcd",
                                                 "dataTakingYear"};
    if (!is_data) {
      output_variables.insert(output_variables.end(), {"sum_weights", "lumi", "xsec", "PileupWeight_NOSYS", "beamSpotWeight", "generatorWeight_NOSYS"});
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

  void Analysis::define_physics_variables()
  {
    const std::string ljet_prefix = "recojet_antikt10UFO_";
    const std::string gn2x_prefix = ljet_prefix + "GN2Xv01_";

    std::vector<std::pair<std::string, std::string>> variables = {
      {"ljet_m", ljet_prefix + "m_NOSYS/1000"},
      {"ljet_pt", ljet_prefix + "pt_NOSYS/1000"},
      {"ljet_phi", ljet_prefix + "phi_NOSYS"},
      {"ljet_eta", ljet_prefix + "eta_NOSYS"},

      {"ljet_phbb", gn2x_prefix + "phbb_NOSYS"},
      {"ljet_phcc", gn2x_prefix + "phcc_NOSYS"},
      {"ljet_ptop", gn2x_prefix + "ptop_NOSYS"},
      {"ljet_pqcd", gn2x_prefix + "pqcd_NOSYS"},

      {"ph_pt", "ph_pt_NOSYS[0]/1000"},
      {"ph_eta", "ph_eta_NOSYS[0]"},
      {"ph_phi", "ph_phi_NOSYS[0]"}
    };

    for (const auto& [name, expr] : variables) {
      df_ = df_.value().Define(name, expr);
    }
  }

  void Analysis::select_Z_candidate()
  {
    // Mass and pt selection
    auto zcand_m_pt_cut = [&](float ljet_m, float ljet_pt) -> bool {
        return ljet_m  > config_.analysis.min_mass &&
               ljet_m  < config_.analysis.max_mass &&
               ljet_pt > config_.analysis.min_pt &&
               ljet_pt < config_.analysis.max_pt;
    };

    // Delta phi selection
    auto zcand_dphi_cut = [&](float ljet_phi, float ph_phi) -> bool {
        float dphi = std::abs(ljet_phi - ph_phi);
        dphi = dphi > M_PI ? 2 * M_PI - dphi : dphi;
        return dphi > config_.analysis.min_dphi;
    };

    //Photon
    df_ = df_.value()
             .Filter("ph_baselineSelection_Tight_FixedCutTight_NOSYS[0]")
             .Filter("ph_passesOR_NOSYS[0]");
    update_cutflow("Photon ISO and ID");

    df_ = df_.value()
             .Filter("ph_pt > 175");

    update_cutflow("Photon pT cut");

    // Apply selection masks and define Z candidate variables
    df_ = df_.value()
             .Define("zcand_m", "ljet_m[0]")
             .Define("zcand_pt", "ljet_pt[0]")
             .Define("zcand_phi", "ljet_phi[0]")
             .Define("zcand_eta", "ljet_eta[0]")
             .Define("zcand_phbb", "ljet_phbb[0]")
             .Define("zcand_phcc", "ljet_phcc[0]")
             .Define("zcand_ptop", "ljet_ptop[0]")
             .Define("zcand_pqcd", "ljet_pqcd[0]");

    df_ = df_.value().Filter(zcand_m_pt_cut, {"zcand_m", "zcand_pt"});
    update_cutflow("Z boson mass and pT selection");
    df_ = df_.value().Filter(zcand_dphi_cut, {"zcand_phi", "ph_phi"});
    update_cutflow("Delta(Z, photon) cut");

    df_ = df_.value()
             .Define("zcand_log_phbb", "-log(zcand_phbb)")
             .Define("zcand_log_phcc", "-log(zcand_phcc)")
             .Define("zcand_log_ptop", "-log(zcand_ptop)")
             .Define("zcand_log_pqcd", "-log(zcand_pqcd)");
  }

  void Analysis::get_dhbb()
  {
    float fcc = 0.02;

    auto make_dhbb_calculator = [fcc](float ftop) {
      return [fcc, ftop](float phbb, float phcc, float ptop, float pqcd) -> float {
        return std::log(phbb/(fcc*phcc + ftop*ptop + (1 - ftop - fcc) * pqcd));
      };
    };

    auto compute_dhbb = make_dhbb_calculator(0.25);
    auto compute_dhbb_ftop0 = make_dhbb_calculator(0.0);

    df_ = df_.value()
             .Define("zcand_gn2x", compute_dhbb,
                    {"zcand_phbb", "zcand_phcc", "zcand_ptop", "zcand_pqcd"})
             .Define("zcand_gn2x_ftop0", compute_dhbb_ftop0,
                    {"zcand_phbb", "zcand_phcc", "zcand_ptop", "zcand_pqcd"});
  }

  void Analysis::save_histograms(const std::string& output_file_name, const std::string& weight, std::string selection)
  {
    std::unique_ptr<TFile> outFile(TFile::Open(output_file_name.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        throw std::runtime_error("Failed to create output file: " + output_file_name);
    }

    spdlog::info("Creating and saving histograms to {}", output_file_name);

    auto df = df_.value().Filter(selection).Filter(weight + " < 1e5");
    auto h_m_zcand = df.Histo1D({"h_m_zcand", "Z candidate mass;m [GeV];Events", 24, 40, 160},
                                  "zcand_m", weight);
    auto h_pt_zcand = df.Histo1D({"h_pt_zcand", "Z candidate p_{T};m [GeV];Events", 25, 200, 500},
                                  "zcand_pt", weight);
    auto h_m_pt_zcand = df.Histo2D({"h_m_pt_zcand", "Z candidate mass and p_{T}", 24, 40, 160, 25, 200, 500},
                                   "zcand_m", "zcand_pt", weight);
    auto h_gn2x = df.Histo1D({"h_gn2x", "DXbb;;Events", 25, -9, 9},
                               "zcand_gn2x", weight);
    auto h_gn2x_ftop0 = df.Histo1D({"h_gn2x_ftop0", "DXbb;;Events", 25, -9, 9},
                                     "zcand_gn2x_ftop0", weight);

    auto h_phbb_zcand = df.Histo1D({"h_phbb_zcand", "-log(phbb);-log(phbb);Events", 20, 0, 0.2},
                                    "zcand_log_phbb", weight);

    auto h_phcc_zcand = df.Histo1D({"h_phcc_zcand", "-log(phcc);-log(phcc);Events", 25, 0, 17},
                                    "zcand_log_phcc", weight);

    auto h_ptop_zcand = df.Histo1D({"h_ptop_zcand", "-log(ptop);-log(ptop);Events", 25, 0, 11},
                                    "zcand_log_ptop", weight);

    auto h_pqcd_zcand = df.Histo1D({"h_pqcd_zcand", "-log(pqcd);-log(pqcd);Events", 25, 0, 8},
                                    "zcand_log_pqcd", weight);

    auto h_m_phbb_zcand = df.Histo2D({"h_m_phbb_zcand", "Z candidate mass and p_{T}", 24, 40, 160, 24, 0, 0.48},
                                   "zcand_m", "zcand_log_phbb", weight);

    auto h_pqcd_phbb_zcand = df.Histo2D({"h_pqcd_phbb_zcand", "Z candidate pqcd and phbb", 80, 0, 8, 90, 0, 9},
                                   "zcand_log_pqcd", "zcand_log_phbb", weight);

    h_m_zcand->Write();
    h_pt_zcand->Write();
    h_m_pt_zcand->Write();
    h_phbb_zcand->Write();
    h_phcc_zcand->Write();
    h_ptop_zcand->Write();
    h_pqcd_zcand->Write();
    h_gn2x->Write();
    h_gn2x_ftop0->Write();
    h_m_phbb_zcand->Write();
    h_pqcd_phbb_zcand->Write();

    for (auto wp : config_.analysis.wps) {
      std::string pass_dhbb = "pass_" + wp;
      std::string histo_name = "h_m_zcand_cut_" + wp;
      std::string histo_name_phbb = "h_phbb_zcand_cut_" + wp;
      std::string histo_name_phcc = "h_phcc_zcand_cut_" + wp;
      std::string histo_name_ptop = "h_ptop_zcand_cut_" + wp;
      std::string histo_name_pqcd = "h_pqcd_zcand_cut_" + wp;
      auto df_cut = df.Filter(pass_dhbb);
      auto h_m_zcand_cut = df_cut.Histo1D({histo_name.c_str(), "Z candidate mass;m [GeV];Events", 24, 40, 160}, "zcand_m", weight);
      auto h_phbb_zcand_cut = df_cut.Histo1D({histo_name_phbb.c_str(), "-log(phbb);-log(phbb);Events", 20, 0, 0.2}, "zcand_log_phbb", weight);
      auto h_phcc_zcand_cut = df_cut.Histo1D({histo_name_phcc.c_str(), "-log(phcc);-log(phcc);Events", 20, 0, 17}, "zcand_log_phcc", weight);
      auto h_ptop_zcand_cut = df_cut.Histo1D({histo_name_ptop.c_str(), "-log(ptop);-log(ptop);Events", 20, 0, 11}, "zcand_log_ptop", weight);
      auto h_pqcd_zcand_cut = df_cut.Histo1D({histo_name_pqcd.c_str(), "-log(pqcd);-log(pqcd);Events", 20, 0, 8}, "zcand_log_pqcd", weight);
      h_m_zcand_cut->Write();
      h_phbb_zcand_cut->Write();
      h_phcc_zcand_cut->Write();
      h_ptop_zcand_cut->Write();
      h_pqcd_zcand_cut->Write();
    }

    int n_cuts = cutflow_.size();
    auto h_cutflow = new TH1D("h_cutflow", "Cutflow Histogram;Cut;Weighted Events", n_cuts, 0, n_cuts);
    auto h_cutflow_unweighted = new TH1D("h_cutflow_unweighted", "Cutflow Histogram;Cut;Unweighted Events", n_cuts, 0, n_cuts);

    for (int i = 0; i < n_cuts; ++i) {
        h_cutflow->SetBinContent(i + 1, cutflow_[i].second);
        h_cutflow_unweighted->SetBinContent(i + 1, cutflow_unweighted_[i].second);
        h_cutflow->GetXaxis()->SetBinLabel(i + 1, cutflow_[i].first.c_str());
        h_cutflow_unweighted->GetXaxis()->SetBinLabel(i + 1, cutflow_[i].first.c_str());
    }
    h_cutflow->Write();
    h_cutflow_unweighted->Write();

    spdlog::info("Successfully wrote histograms to output file");
  }
}
