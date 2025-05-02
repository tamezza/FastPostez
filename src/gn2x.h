#pragma once

#include <vector>
#include <string>
#include <map>



class GN2XHandler {
  struct GN2XCut {
    double min_mass;
    double max_mass;
    double cut;
  };

  private:
    std::map<std::string, std::vector<GN2XCut>> all_cuts_;
  public:
    GN2XHandler(const std::string wps_cuts_filename);
    std::string make_selection_code(const std::string& wp_name,
                                    const std::string& gn2x_var,
                                    const std::string& mass_var);
    std::vector<std::string> get_wps();

  private:
    void load_flatmass_cuts_CSV(const std::string& filename);
};
