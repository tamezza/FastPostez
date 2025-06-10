#pragma once

#include <vector>
#include <string>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <nlohmann/json.hpp>


class GN2XHandler {
  struct FlatMass {
    std::vector<int> mass;
    std::vector<float> cutvalues;
  };
  typedef std::pair<std::vector<int>, std::vector<FlatMass>> MassPtCuts;

  private:
    std::map<std::string, MassPtCuts> all_cuts_;
  public:
    explicit GN2XHandler(const std::string& wps_cuts_filename);
    ROOT::RDF::RNode define_pass(ROOT::RDF::RNode df,
                                 const std::string& dhbb,
                                 const std::string& ljet_m,
                                 const std::string& ljet_pt);

  private:
    void load_flatmass_cuts(const std::string& xbb_file_name,
                            const std::string& tagger,
                            const std::string& large_jet);
    MassPtCuts get_mass_pt_cuts(nlohmann::json j);
    FlatMass get_flat_mass(nlohmann::json j);
    float get_cut(const MassPtCuts& cuts, float mass, float pt);
    int find_pt_bin(const std::vector<int>& bin_edges, float pt);
    int find_mass_bin(const std::vector<int>& bin_edges, float mass);
};
