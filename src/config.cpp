#include <string>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include "config.h"

config::Config config::load_config(const std::string& config_file_name)
{
  config::Config config;

  std::ifstream config_file(config_file_name);
  if (!config_file.is_open()) {
    throw std::runtime_error("Could not open config file: " + config_file_name);
  }

  nlohmann::json j;
  try {
      config_file >> j;
  } catch (const nlohmann::json::parse_error& e) {
      throw std::runtime_error("Failed to parse config file '" + config_file_name + "': " + e.what());
  }

  try {
    config.ntuple.tree_name  = j.at("ntuple").at("tree_name").get<std::string>();

    config.analysis.dsids    = j.at("analysis").at("dsids").get<std::vector<int>>();
    config.analysis.years    = j.at("analysis").at("years").get<std::vector<int>>();
    config.analysis.metadata = j.at("analysis").at("metadata").get<std::string>();
    config.analysis.flatmass = j.at("analysis").at("flatmass").get<std::string>();
    config.analysis.wps      = j.at("analysis").at("wps").get<std::vector<std::string>>();

    config.analysis.min_mass = j.at("analysis").at("min_mass").get<float>();
    config.analysis.max_mass = j.at("analysis").at("max_mass").get<float>();
    config.analysis.min_pt   = j.at("analysis").at("min_pt").get<float>();
    config.analysis.max_pt   = j.at("analysis").at("max_pt").get<float>();
    config.analysis.min_dphi = j.at("analysis").at("min_dphi").get<float>();

  } catch (const nlohmann::json::out_of_range& e) {
      throw std::runtime_error(std::string("Missing config key in JSON: ") + e.what());
  } catch (const nlohmann::json::type_error& e) {
      throw std::runtime_error(std::string("Config type mismatch in JSON: ") + e.what());
  }

  return config;
}
