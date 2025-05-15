#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TChain.h>
#include <spdlog/spdlog.h>
#include "output_manager.h"
#include "input_handler.h"
#include "weights.h"
#include "config.h"
#include "xbby_calib.h"
#include "gn2x.h"

namespace xbbycalib {
  Analysis::Analysis(const std::string config_file,
                     const std::string input_folder,
                     std::string output_folder,
                     bool multi_threading)
    : input_folder_(input_folder),
      output_folder_(output_folder),
      multi_threading_(multi_threading)
  {
    spdlog::info("Initializing analysis...");
    spdlog::info("Config file: {}", config_file);
    spdlog::info("Input folder: {}", input_folder);
    spdlog::info("Output folder: {}", output_folder);

    spdlog::info("Loading config from {}", config_file);
    config_ = config::load_config(config_file);
  }

  void Analysis::run()
  {
    output_folder_ = output_manager::create_output_folder(output_folder_);
    output_manager::set_logger(output_folder_);
    if (multi_threading_) {
      ROOT::EnableImplicitMT();
      spdlog::info("Multi-threading enabled.");
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

    define_physics_variables();
    select_Z_candidate();
    get_dhbb();
    GN2XHandler gn2x_handler_(config_.analysis.flatmass);

    for (auto wp : config_.analysis.wps) {
      auto dhbb_selection_code = gn2x_handler_.make_selection_code(wp, "zcand_gn2x", "zcand_m");
      std::string pass_dhbb = "pass_dhbb_" + wp;
      df_ = df_.value()
               .Define(pass_dhbb, dhbb_selection_code)
               .Define("zcand_m_" + pass_dhbb, "zcand_m * " + pass_dhbb);
    }
    const std::string output_file_name = output_folder_ + "/histograms/hists_"
                                                        + std::to_string(sample_label) + ".root";
    save_histograms(output_file_name);
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
    auto zcand_m_pt_mask = [&](const ROOT::RVecF& ljet_m, const ROOT::RVecF& ljet_pt) -> ROOT::RVecB {
        return ljet_m  > config_.analysis.min_mass &&
               ljet_m  < config_.analysis.max_mass &&
               ljet_pt > config_.analysis.min_pt &&
               ljet_pt < config_.analysis.max_pt;
    };

    // Delta phi selection
    auto zcand_dphi_mask = [&](const ROOT::RVecF& ljet_phi, float ph_phi) -> ROOT::RVecB {
        return ROOT::VecOps::Map(ljet_phi, [ph_phi](float phi) {
            float dphi = std::abs(phi - ph_phi);
            return dphi > M_PI ? 2 * M_PI - dphi : dphi;
        }) > config_.analysis.min_dphi;
    };

    // Apply selection masks and define Z candidate variables
    df_ = df_.value()
             .Define("zcand_m_pt_mask", zcand_m_pt_mask, {"ljet_m", "ljet_pt"})
             .Define("zcand_dphi_mask", zcand_dphi_mask, {"ljet_phi", "ph_phi"})
             .Define("zcand_mask", "zcand_m_pt_mask && zcand_dphi_mask")
             .Define("zcand_m", "ljet_m[zcand_mask][0]")
             .Define("zcand_pt", "ljet_pt[zcand_mask][0]")
             .Define("zcand_phi", "ljet_phi[zcand_mask][0]")
             .Define("zcand_eta", "ljet_eta[zcand_mask][0]")
             .Define("zcand_phbb", "ljet_phbb[zcand_mask][0]")
             .Define("zcand_phcc", "ljet_phcc[zcand_mask][0]")
             .Define("zcand_ptop", "ljet_ptop[zcand_mask][0]")
             .Define("zcand_pqcd", "ljet_pqcd[zcand_mask][0]")
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

  void Analysis::save_histograms(const std::string& output_file_name)
  {
    std::unique_ptr<TFile> outFile(TFile::Open(output_file_name.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        throw std::runtime_error("Failed to create output file: " + output_file_name);
    }

    spdlog::info("Creating and saving histograms to {}", output_file_name);

    auto& df = df_.value();
    auto h_m_zcand = df.Histo1D({"h_m_zcand", "Z candidate mass;m [GeV];Events", 24, 40, 160},
                                  "zcand_m", "total_weight");
    auto h_gn2x = df.Histo1D({"h_m_gn2x", "DXbb;;Events", 20, -10, 10},
                               "zcand_gn2x", "total_weight");
    auto h_gn2x_ftop0 = df.Histo1D({"h_m_gn2x_ftop0", "DXbb;;Events", 20, -10, 10},
                                     "zcand_gn2x_ftop0", "total_weight");

    auto h_phbb_zcand = df.Histo1D({"h_phbb_zcand", "log(phbb);log(phbb);Events", 50, 0, 10},
                                    "zcand_log_phbb", "total_weight");

    auto h_phcc_zcand = df.Histo1D({"h_phcc_zcand", "log(phcc);log(phcc);Events", 50, 0, 20},
                                    "zcand_log_phcc", "total_weight");

    auto h_ptop_zcand = df.Histo1D({"h_ptop_zcand", "log(ptop);log(ptop);Events", 50, 0, 12},
                                    "zcand_log_ptop", "total_weight");

    auto h_pqcd_zcand = df.Histo1D({"h_pqcd_zcand", "log(pqcd);log(pqcd);Events", 50, 0, 10},
                                    "zcand_log_pqcd", "total_weight");

    auto h_phbb_zcand_sidebands = df.Filter("zcand_m < 65 || zcand_m > 110").Histo1D({"h_phbb_zcand_sidebands", "log(phbb);log(phbb);Events", 50, 0, 10},
                                    "zcand_log_phbb", "total_weight");

    auto h_phcc_zcand_sidebands = df.Filter("zcand_m < 65 || zcand_m > 110").Histo1D({"h_phcc_zcand_sidebands", "log(phcc);log(phcc);Events", 50, 0, 20},
                                    "zcand_log_phcc", "total_weight");

    auto h_ptop_zcand_sidebands = df.Filter("zcand_m < 65 || zcand_m > 110").Histo1D({"h_ptop_zcand_sidebands", "log(ptop);log(ptop);Events", 50, 0, 12},
                                    "zcand_log_ptop", "total_weight");

    auto h_pqcd_zcand_sidebands = df.Filter("zcand_m < 65 || zcand_m > 110").Histo1D({"h_pqcd_zcand_sidebands", "log(pqcd);log(pqcd);Events", 50, 0, 10},
                                    "zcand_log_pqcd", "total_weight");

    h_m_zcand->Write();
    h_phbb_zcand->Write();
    h_phcc_zcand->Write();
    h_ptop_zcand->Write();
    h_pqcd_zcand->Write();
    h_phbb_zcand_sidebands->Write();
    h_phcc_zcand_sidebands->Write();
    h_ptop_zcand_sidebands->Write();
    h_pqcd_zcand_sidebands->Write();
    h_gn2x->Write();
    h_gn2x_ftop0->Write();

    for (auto wp : config_.analysis.wps) {
      std::string zcand_pass_dhbb = "zcand_m_pass_dhbb_" + wp;
      std::string histo_name = "h_m_zcand_cut_" + wp;
      auto h_m_zcand_cut = df.Histo1D({histo_name.c_str(), "Z candidate mass;m [GeV];Events", 24, 40, 160},
                                        zcand_pass_dhbb.c_str(), "total_weight");
      h_m_zcand_cut->Write();
    }

    spdlog::info("Successfully wrote histograms to output file");
  }
}
