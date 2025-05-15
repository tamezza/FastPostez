#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <map>
#include <algorithm>
#include <spdlog/spdlog.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObject.h>
#include <TCollection.h>
#include <ROOT/RDataFrame.hxx>

#include "weights.h"
#include "input_handler.h"

namespace fs = std::filesystem;

InputHandler::InputHandler(const std::string& input_folder,
                           const std::string& tree_name,
                           int sample_label,
                           bool is_data,
                           const std::string& metadata)
  : input_folder_(input_folder),
    tree_name_(tree_name),
    sample_label_(sample_label),
    is_data_(is_data),
    metadata_(metadata)
{
}

std::unique_ptr<TChain> InputHandler::get_chain()
{
  std::string label = std::to_string(sample_label_);
  std::string pattern = is_data_ ? "grp" + label : label;
  const auto input_folders_files = InputHandler::make_file_list(pattern);

  all_input_files_.clear();
  for (const auto& [folder, files] : input_folders_files) {
    for (const auto& file : input_folders_files.at(folder)) {
      if (file.find(merged_file_name_) == std::string::npos) {
        all_input_files_.emplace_back(file);
      }
    }
  }

  const auto merged_files = make_merged_file_list(input_folders_files);
  return create_safe_chain(merged_files);
}

std::map<std::string, std::vector<std::string>> InputHandler::make_file_list(const std::string& pattern)
{
  std::map<std::string, std::vector<std::string>> folder_file_map;
  std::regex folder_pattern("user.*" + pattern + ".*");

  try {
    for (const auto& folder : fs::directory_iterator(input_folder_)) {
      const std::string folder_name = folder.path().filename().string();
      if(std::regex_search(folder_name, folder_pattern)) {
        std::vector<std::string> file_list;
        for (const auto& file : fs::directory_iterator(folder.path())) {
          file_list.push_back(file.path().string());
        }
        folder_file_map[folder_name] = std::move(file_list);
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return folder_file_map;
}

std::vector<std::string> InputHandler::make_merged_file_list(std::map<std::string, std::vector<std::string>> input_folders_files)
{
  std::vector<std::string> input_files;
  for (const auto& [folder, files] : input_folders_files) {
    input_files.emplace_back(input_folder_ + "/" + folder + "/" + merged_file_name_);
    auto it = std::find_if(files.begin(), files.end(), [&](const std::string& f) {
      return f.find(merged_file_name_) != std::string::npos;
    });
    if (it == files.end()) {
      InputHandler::create_merged_file(folder, files);
    }
  }
  return input_files;
}

void InputHandler::create_merged_file(const std::string& folder, const std::vector<std::string>& files)
{
  auto chain_ptr = InputHandler::create_safe_chain(files);
  ROOT::RDataFrame df(*chain_ptr);
  spdlog::info("Merging input files...");
  auto df_weighted = weights::define_event_weights(df, metadata_, sample_label_, all_input_files_, is_data_);
  df_weighted.Snapshot(tree_name_, input_folder_ + "/" + folder + "/" + merged_file_name_);
}

std::unique_ptr<TChain> InputHandler::create_safe_chain(const std::vector<std::string>& files)
{
  auto chain = std::make_unique<TChain>(tree_name_.c_str());
  std::vector<std::string> valid_files;

  spdlog::info("Processing {} files...", files.size());

  int total_files = files.size();
  int processed_files = 0;

  for (const auto& file_path : files) {
    std::unique_ptr<TFile> file(TFile::Open(file_path.c_str()));
    processed_files++;

    if (!file || file->IsZombie()) {
      spdlog::warn("Could not open {}", file_path);
      continue;
    }

    auto tree = dynamic_cast<TTree*>(file->Get(tree_name_.c_str()));
    if (!tree) {
      spdlog::warn("No tree '{}' in {}", tree_name_, file_path);
      continue;
    }

    if (!tree->GetEntries()) continue;

    valid_files.push_back(file_path);
    chain->Add(file_path.c_str());

    float progress = static_cast<float>(processed_files) / total_files;
    std::cout << "Building the list of valid input files... "
              << std::fixed << std::setprecision(0)
              << progress * 100.0 << "%\r";
    std::cout.flush();
  }
  std::cout.flush();
  spdlog::info("Successfully added {} valid files to the chain.", valid_files.size());

  if (!chain || chain->GetEntries() <= 0) {
    throw std::runtime_error("Failed to create valid TChain from input files");
  }

  spdlog::info("Created chain with {} entries", chain->GetEntries());

  return chain;
}

