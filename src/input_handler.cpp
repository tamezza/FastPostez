#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <map>
#include <spdlog/spdlog.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObject.h>
#include <TCollection.h>
#include <ROOT/RVec.hxx>
#include "weights.h"
#include "input_handler.h"

namespace fs = std::filesystem;

std::vector<std::string> input_handler::make_file_list(const std::string& input_folder, const std::string& pattern)
{
  std::vector<std::string> file_list;
  std::regex folder_pattern("user.*" + pattern + ".*");

  try {
    for (const auto& folder : fs::directory_iterator(input_folder)) {
      if (folder.is_directory() && std::regex_search(folder.path().filename().string(), folder_pattern)) {
        std::transform(fs::directory_iterator(folder.path()),
            fs::directory_iterator(),
            std::back_inserter(file_list),
            [](const auto& file) { return file.path().string(); });
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return file_list;
}

std::vector<std::string> input_handler::prepare_input_files(const std::string& input_folder,
                                                            int sample_label, bool is_data)
{
  spdlog::info("Building list of input files...");
  std::string label = std::to_string(sample_label);
  std::string pattern = is_data ? "grp" + label : label;
  std::vector<std::string> input_files = input_handler::make_file_list(input_folder, pattern);

  if (input_files.empty()) {
        throw std::runtime_error("No input files found for dataset ID: " + label);
  }

  spdlog::info("Found {} input files", input_files.size());
  return input_files;
}

std::unique_ptr<TChain> input_handler::create_safe_chain(const std::vector<std::string>& file_list,
                                                         const std::string& tree_name)
{
  auto chain = std::make_unique<TChain>(tree_name.c_str());
  std::vector<std::string> valid_files;

  spdlog::info("Processing {} files...", file_list.size());

  int total_files = file_list.size();
  int processed_files = 0;

  for (const auto& file_path : file_list) {
    std::unique_ptr<TFile> file(TFile::Open(file_path.c_str()));
    processed_files++;

    if (!file || file->IsZombie()) {
      spdlog::warn("Could not open {}", file_path);
      continue;
    }

    auto tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
    if (!tree) {
      spdlog::warn("No tree '{}' in {}", tree_name, file_path);
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

void merge_input_files(const std::vector<std::string>& file_list, const std::string& tree_name, int sample_label)
{
  fs::path output_path = output_path_str;
}
