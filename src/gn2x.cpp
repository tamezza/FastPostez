#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <spdlog/spdlog.h>
#include "gn2x.h"

GN2XHandler::GN2XHandler(const std::string wps_cuts_filename)
{
  load_flatmass_cuts_CSV(wps_cuts_filename);
}

void GN2XHandler::load_flatmass_cuts_CSV(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open()) {
      throw std::runtime_error("Error: Could not open CSV file: " + filename);
    }

    std::string line;
    // Store the mapping of column index to working point name
    std::vector<std::pair<int, std::string>> wp_columns;

    // --- Read Header ---
    if (std::getline(file, line)) {
      std::stringstream header_stream(line);
      std::string cell;
      int col_index = 0;
      while (std::getline(header_stream, cell, ',')) {
        // Trim whitespace
        cell.erase(0, cell.find_first_not_of(" \t\n\r\f\v"));
        cell.erase(cell.find_last_not_of(" \t\n\r\f\v") + 1);

        wp_columns.push_back({col_index, cell});
        all_cuts_[cell] = {}; // Initialize vector for this WP

        col_index++;
      }
    } else {
      throw std::runtime_error("Error: CSV file is empty or could not read header: " + filename);
    }

    // --- Read Data Rows ---
    int line_number = 1; // Start counting after header
    while (std::getline(file, line)) {
      line_number++;
      // Skip empty lines or lines starting with # (comments)
      if (line.empty() || line[0] == '#') continue;

      std::stringstream line_stream(line);
      std::string cell;
      std::vector<std::string> cells;

      while (std::getline(line_stream, cell, ',')) {
        cells.push_back(cell);
      }

      // --- Parse Mass Range  ---
      double min_mass = 0;
      double max_mass = 0;
      std::sscanf(cells[0].c_str(), "( M()>= %lf && M() < %lf )", &min_mass, &max_mass);

      // --- Parse Cut Values for each Working Point ---
      for (const auto& wp_info : wp_columns) {
        int wp_col_index = wp_info.first;
        const std::string& wp_name = wp_info.second;
        if (wp_name == "WP") continue;

        std::string value_str = cells[wp_col_index];
        //spdlog::info("wp: {}", wp_name);
        //spdlog::info("value: {}", value_str);
        double cut_value = std::stod(value_str);
        all_cuts_[wp_name].push_back({min_mass, max_mass, cut_value});
      }
    }
    file.close();
}

std::string GN2XHandler::make_selection_code(const std::string& wp_name,
                                             const std::string& gn2x_var,
                                             const std::string& mass_var)
{
  auto cuts = all_cuts_[wp_name];
  std::string selection_code = "if (" + mass_var + " < 40) return false;\n";
  for (const auto cut : cuts) {
    selection_code += "else if (" + mass_var + " >= " + std::to_string(cut.min_mass)
                    + " && "      + mass_var + " < "  + std::to_string(cut.max_mass)
                    + ") {return " + gn2x_var + " > " + std::to_string(cut.cut)
                    + ";}\n";
  }
  selection_code += "else return false;\n";
  return selection_code;
}

std::vector<std::string> GN2XHandler::get_wps()
{
  std::vector<std::string> wps;
  for (const auto wp : all_cuts_) {
    wps.push_back(wp.first);
  }
  return wps;
}
