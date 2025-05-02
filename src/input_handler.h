#pragma once

#include <string>
#include <vector>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <ROOT/RVec.hxx>

namespace input_handler {
  std::vector<std::string> make_file_list(const std::string& input_folder, const std::string& pattern);
  std::vector<std::string> prepare_input_files(const std::string& input_folder, int sample_label, bool is_data);
  std::unique_ptr<TChain> create_safe_chain(const std::vector<std::string>& file_list, const std::string& tree_name);
}
