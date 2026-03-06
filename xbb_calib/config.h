#pragma once

#include <string>
#include <vector>
#include <map>

namespace config {
  struct Config {
    struct {
      std::string tree_name;
    } ntuple;
    struct {
      std::vector<int> dsids;
      std::vector<int> years;
      std::map<int, std::vector<std::string>> trigger_map;
      std::string metadata;
      std::string flatmass;
      std::vector<std::string> wps;

      float min_mass;
      float max_mass;
      float min_pt;
      float max_pt;
      float min_dphi;

    } analysis;
  };

  config::Config load_config(const std::string& config_file_name);
}

