#pragma once

#include <string>
#include <vector>
#include <optional>
#include <ROOT/RDataFrame.hxx>
#include "config.h"

namespace xbbycalib {
  class Analysis {
    private:
      const std::string input_folder_;
      std::string output_folder_;
      bool multi_threading_;
      config::Config config_;
      std::optional<ROOT::RDF::RNode> df_;

    public:
      Analysis(const std::string config_file,
               const std::string input_folder,
               std::string output_folder,
               bool multi_threading);
      void run();

    private:
      void run_sample(int sample_label, bool is_data);
      void define_physics_variables();
      void select_Z_candidate();
      void get_dhbb();
      std::string create_dhbb_selection_code(const std::string& csv_file,
                                             const std::string& wp_name,
                                             const std::string& gn2x_var);
      void save_histograms(const std::string& output_file_name);
  };
}
