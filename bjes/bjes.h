#pragma once

#include <string>
#include <vector>
#include <optional>
#include <ROOT/RDataFrame.hxx>
#include "config_bjes.h"

namespace bjes {
  class Analysis {
    private:
      const std::string input_folder_;
      std::string output_folder_;
      int n_threads_;
      config_bjes::Config config_;
      std::optional<ROOT::RDF::RNode> df_;
      std::vector<std::pair<std::string, double>> cutflow_;
      std::vector<std::pair<std::string, double>> cutflow_unweighted_;

    public:
      Analysis(const config_bjes::Config& config,
               const std::string& input_folder,
               std::string output_folder,
               int n_threads);
      void run();

    private:
      void run_sample(int sample_label, bool is_data);
      void update_cutflow(std::string label);
      void apply_triggers();
      void define_physics_variables(bool is_data);
      void get_dhbb();
      std::string create_dhbb_selection_code(const std::string& csv_file,
                                             const std::string& wp_name,
                                             const std::string& gn2x_var);
      std::string get_branch_name(const std::string& base_name,
                                  const std::vector<std::string>& candidates);
  };
}
