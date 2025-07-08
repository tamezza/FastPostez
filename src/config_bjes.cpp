#include <string>
#include <map>
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include "config_bjes.h"

config_bjes::Config config_bjes::load_config(const std::string& config_file_name)
{
  Config config;

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
    config.ntuple.tree_name = j.at("ntuple").at("tree_name").get<std::string>();

    config.analysis.dsids       = j.at("analysis").at("dsids").get<std::vector<int>>();
    config.analysis.years       = j.at("analysis").at("years").get<std::vector<int>>();
    for (const auto& [year_str, triggers] : j.at("analysis").at("trigger_map").items()) {
      int year = std::stoi(year_str);
      config.analysis.trigger_map[year] = triggers.get<std::vector<std::string>>();
    }
    config.analysis.metadata    = j.at("analysis").at("metadata").get<std::string>();
    config.analysis.analysis    = j.at("analysis").at("analysis").get<std::string>();

  } catch (const nlohmann::json::out_of_range& e) {
      throw std::runtime_error(std::string("Missing config key in JSON: ") + e.what());
  } catch (const nlohmann::json::type_error& e) {
      throw std::runtime_error(std::string("Config type mismatch in JSON: ") + e.what());
  }

  return config;
}
