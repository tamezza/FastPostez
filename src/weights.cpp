#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <regex>
#include <map>
#include <cmath>
#include <spdlog/spdlog.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObject.h>
#include <TCollection.h>
#include <ROOT/RVec.hxx>
#include "weights.h"

namespace weights {
  std::map<int, double> load_metadata(const std::string& filename)
  {
    std::map<int, double> map_sx;
    TTree metadata;

    metadata.ReadFile(filename.c_str());

    int dataset_number;
    double crossSection_pb, kFactor, genFiltEff;
    metadata.SetBranchAddress("dataset_number", &dataset_number);
    metadata.SetBranchAddress("crossSection_pb", &crossSection_pb);
    metadata.SetBranchAddress("kFactor", &kFactor);
    metadata.SetBranchAddress("genFiltEff", &genFiltEff);

    for (int i = 0; i < metadata.GetEntries(); ++i) {
      metadata.GetEntry(i);
      map_sx[dataset_number] = crossSection_pb * kFactor * genFiltEff;
    }

    return map_sx;
  }

  double get_xcross(int dsid, const std::string& metadata_file)
  {
    static std::map<int, double> map_sx = load_metadata(metadata_file);

    auto it = map_sx.find(dsid);
    return (it != map_sx.end()) ? it->second : -1.0;
  }

  double get_sum_w(const std::vector<std::string>& input_files)
  {
    TH1F* merged_hist = nullptr;
    std::regex pattern("^CutBookkeeper_.*_.*_NOSYS$");

    for (const auto& file_path : input_files) {
      auto root_file = TFile::Open(file_path.c_str(), "READ");
      if (!root_file || root_file->IsZombie()) {
        std::cerr << "Warning: Could not open " << file_path << std::endl;
        delete root_file;
        continue;
      }

      TIter next(root_file->GetListOfKeys());
      TKey* key;
      while (key = static_cast<TKey*>(next())) {
        TObject* obj = key->ReadObj();
        auto hist = dynamic_cast<TH1F*>(obj);
        if (hist && std::regex_match(hist->GetName(), pattern)) {
          if (!merged_hist) {
            merged_hist = static_cast<TH1F*>(hist->Clone("MergedHistogram"));
            merged_hist->SetDirectory(nullptr); // Detach from file
          } else {
            merged_hist->Add(hist);
          }
        }
      }

      root_file->Close();
      delete root_file;
    }

    return merged_hist->GetBinContent(2);
  }

  std::vector<std::string> checkWeights(ROOT::RDF::RNode df)
  {
    std::vector<std::string> colNames = df.GetColumnNames();
    std::vector<std::string> listWeights;
    std::vector<std::string> weights = {"PileupWeight_NOSYS", "beamSpotWeight",
      "weight_leptonSF", "weight_photonSF"};

    for (const auto& weight : weights) {
      std::string weightExpr = (std::find(colNames.begin(), colNames.end(), weight) != colNames.end())
        ? weight : "1";

      if (weightExpr == "1") spdlog::warn("Weight: {} not found, setting it to one.", weight);

      listWeights.push_back(weightExpr);
    }
    return listWeights;
  }

  std::string getWeightExpr(double xcross, double sumWeights, std::string period, ROOT::RDF::RNode df)
  {
    // https://twiki.cern.ch/twiki/bin/view/Atlas/LuminosityForPhysics#2015_2018_13_TeV_proton_proton_f
    const std::unordered_map<std::string, double> lumiMap = {
      {"run2", 140068.94},
      {"r16167", 3244.54 + 33402.2},
      {"r13144", 44630.6},
      {"r14145", 58791.6}
    };

    auto it = lumiMap.find(period);
    double lumi = (it != lumiMap.end()) ? it->second : 0.0;

    auto weights = checkWeights(df);
    std::string weightPileup = weights[0];
    std::string weightBeamspot = weights[1];
    std::string weightLeptonSF = weights[2];
    std::string weightPhotonSF = weights[3];

    std::string total_weight = "generatorWeight_NOSYS * " + weightPileup + " * " + weightBeamspot + " * " +
      weightLeptonSF + " * " + weightPhotonSF + " * " +
      std::to_string(lumi) + " * " + std::to_string(xcross) + " / " +
      std::to_string(sumWeights);
    spdlog::info("Weight expression: {}", total_weight);
    return total_weight;
  }

  ROOT::RDF::RNode define_event_weights(ROOT::RDF::RNode df, const std::string& metadata, int sample_label,
                                        const std::vector<std::string>& input_files, bool is_data)
  {
    if (is_data) return df.Define("total_weight", "1");

    const double cross_section = weights::get_xcross(sample_label, metadata);
    const double sum_weights = weights::get_sum_w(input_files);
    const std::string weight_expr = weights::getWeightExpr(cross_section, sum_weights, "run2", df);

    return df.Define("total_weight", weight_expr);
  }
}
