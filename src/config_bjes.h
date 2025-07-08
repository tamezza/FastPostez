#pragma once

#include <string>
#include <vector>
#include <map>

namespace config_bjes {
  struct Config {
    struct {
      std::string tree_name;
    } ntuple;
    struct {
      std::vector<int> dsids;
      std::vector<int> years;
      std::map<int, std::vector<std::string>> trigger_map;
      std::string metadata;
      std::string analysis;
    } analysis;
  };

  Config load_config(const std::string& config_file_name);
}

