#pragma once

#include <string>
#include <vector>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <ROOT/RVec.hxx>

namespace weights {
  std::map<int, double> load_metadata(const std::string& filename);
  double get_xcross(int dsid, const std::string& metadata_file);
  double get_sum_w(const std::vector<std::string>& input_files);
  std::vector<std::string> checkWeights(ROOT::RDF::RNode df);
  std::string getWeightExpr(double xcross, double sumWeights, std::string period, ROOT::RDF::RNode df);
  ROOT::RDF::RNode define_event_weights(ROOT::RDF::RNode df, const std::string& metadata, int sample_label,
                                        const std::vector<std::string>& input_files, bool is_data);
}
