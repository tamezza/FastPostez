#pragma once

#include <string>
#include <vector>
#include <map>
#include <TChain.h>

class InputHandler {
  private:
    const std::string input_folder_;
    const std::string tree_name_;
    int sample_label_;
    bool is_data_;
    const std::string metadata_;
    std::vector<std::string> all_input_files_;
    const std::string merged_file_name_ = "merged.root";
  public:
    InputHandler(const std::string& input_folder, const std::string& tree_name, int sample_label, bool is_data, const std::string& metadata);
    std::unique_ptr<TChain> get_chain();
  private:
    std::map<std::string, std::vector<std::string>> make_file_list(const std::string& pattern);
    std::vector<std::string> make_merged_file_list(std::map<std::string, std::vector<std::string>> input_folders_files);
    void create_merged_file(const std::string& folder, const std::vector<std::string>& files);
    std::unique_ptr<TChain> create_safe_chain(const std::vector<std::string>& file_list);
};
