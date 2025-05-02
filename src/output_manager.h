#pragma once

#include <string>

namespace output_manager {
  std::string create_output_folder(const std::string& output_folder);
  std::string get_timestamp();
  void set_logger(const std::string& output_folder);
}
