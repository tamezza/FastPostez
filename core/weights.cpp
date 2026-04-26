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

  std::map<unsigned int, double> get_sum_w(const std::vector<std::string>& input_files)
  {
    std::regex pattern(R"(^CutBookkeeper_.*_(\d+)_NOSYS$)");
    std::map<unsigned int, double> year_sum;

    for (const auto& file_path : input_files) {
      auto root_file = TFile::Open(file_path.c_str(), "READ");
      if (!root_file || root_file->IsZombie()) {
        std::cerr << "Warning: Could not open " << file_path << std::endl;
        delete root_file;
        continue;
      }

      TIter next(root_file->GetListOfKeys());
      TKey* key;
      while ((key = static_cast<TKey*>(next()))) {
        TObject* obj = key->ReadObj();
        auto hist = dynamic_cast<TH1F*>(obj);
        if (!hist) continue;

        const std::string name(hist->GetName());
        std::smatch match;
        if (std::regex_match(name, match, pattern) && match.size() > 1) {
          int run_number = std::stoi(match[1].str());
          double content = hist->GetBinContent(2);

          // Run-3 (13.6 TeV) - MC23 samples
          if (run_number == 470000) {
            year_sum[2024] += content;
          }
          else if (run_number == 450000) {
            year_sum[2023] += content;
          }
          else if (run_number == 410000) {
            year_sum[2022] += content;
          }
          // // Run-2 (13 TeV) - MC20 samples
          else if (run_number == 310000) {
            year_sum[2018] += content;
          }
          else if (run_number == 300000) {
            year_sum[2017] += content;
          }
          else if (run_number == 284500) {
            year_sum[2016] += content;
            year_sum[2015] += content;
          }
          else {
            std::cerr << "Warning: Unrecognized run number " << run_number << " in " << name << std::endl;
            continue;
          }
        }
      }

      root_file->Close();
      delete root_file;
    }

    return year_sum;
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
  
  // Integrated Luminosity pb^(-1)
  double get_lumi(unsigned int year)
  {
    // https://twiki.cern.ch/twiki/bin/view/Atlas/LuminosityForPhysics#2015_2018_13_TeV_proton_proton_f
    const std::unordered_map<int, double> lumiMap = {
      // Run-2: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/GoodRunListsForAnalysisRun2
      //{2015, 3244.54 + 33402.2},
      //{2016, 3244.54 + 33402.2},
      {2015, 3244.54},
      {2016, 33402.2},
      {2017, 44630.6},
      {2018, 58791.6},
      // Run-3: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/GoodRunListsForAnalysisRun3
      {2022, 26328.8}, //pb-1
      {2023, 25204.3}, //pb-1
      {2024, 107900.0} //I mistake here I took the fb-1; {2024, 107.9} or 109.4 fb-1
      
    };

    auto it = lumiMap.find(year);
    return (it != lumiMap.end()) ? it->second : 0.0;
  }

  std::string getWeightExpr(ROOT::RDF::RNode df)
  {
    auto weights = checkWeights(df);
    std::string weightPileup = weights[0];
    //std::string weightBeamspot = weights[1];
    std::string weightLeptonSF = weights[2];
    std::string weightPhotonSF = weights[3];

    //std::string total_weight = "generatorWeight_NOSYS * " + weightPileup + " * " + weightBeamspot + " * " +
    //  weightLeptonSF + " * " + weightPhotonSF + " * xsec * lumi / sum_weights";
    std::string total_weight = "generatorWeight_NOSYS * " + weightPileup + " * " +
      weightLeptonSF + " * " + weightPhotonSF + " * xsec * lumi / sum_weights";
    spdlog::info("Weight expression: {}", total_weight);
    return total_weight;
  }

  ROOT::RDF::RNode define_event_weights(ROOT::RDF::RNode df, const std::string& metadata, int sample_label,
                                        const std::vector<std::string>& input_files, bool is_data)
  {
    if (is_data) return df.Define("total_weight", "1");

    const double cross_section = get_xcross(sample_label, metadata);
    const auto sum_weights = get_sum_w(input_files);

    auto sum_weights_year = [sum_weights](unsigned int year) {return sum_weights.at(year);};

    df = df.Define("xsec", std::to_string(cross_section))
           .Define("sum_weights", sum_weights_year, {"dataTakingYear"})
           .Define("lumi", get_lumi, {"dataTakingYear"});

    return df.Define("total_weight", getWeightExpr(df));
  }
}
