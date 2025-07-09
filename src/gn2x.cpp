#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <ranges>
#include <ROOT/RDataFrame.hxx>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include "gn2x.h"

GN2XHandler::GN2XHandler(const std::string& wps_cuts_filename)
{
  const std::string tagger = "GN2Xv01";
  const std::string large_jet = "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets";

  spdlog::info("Loading flat mass cuts for {}.", tagger);
  load_flatmass_cuts(wps_cuts_filename, tagger, large_jet);
}

void GN2XHandler::load_flatmass_cuts(const std::string& xbb_file_name,
                                     const std::string& tagger,
                                     const std::string& large_jet)
{
  std::ifstream xbb_file(xbb_file_name);

  nlohmann::json j;
  xbb_file >> j;

  const auto content = j.at(tagger).at(large_jet);
  const auto wps = content.at("meta").at("OperatingPoints").get<std::vector<std::string>>();

  for (const auto& wp : wps) {
    all_cuts_[wp] = get_mass_pt_cuts(content.at(wp));
  }
}

GN2XHandler::MassPtCuts GN2XHandler::get_mass_pt_cuts(nlohmann::json j)
{
  const auto content = j.at("pT_mass_2d_cutvalue");
  std::vector<int> pt_bins;
  const auto arr = content.at("pTbins");

  pt_bins.reserve(arr.size());
  for (size_t i = 0; i < arr.size() - 1; ++i) {
    pt_bins.push_back(arr[i].get<int>());
  }

  std::vector<FlatMass> flatmass_cuts;
  for (size_t i = 0; i < pt_bins.size(); ++i) {
    std::string pt_range;
    if (i == pt_bins.size() - 1) {
      pt_range = "pT_" + std::to_string(pt_bins.back()) + "_inf";
    }
    else {
      pt_range = "pT_" + std::to_string(pt_bins[i]) + "_" + std::to_string(pt_bins[i + 1]);
    }

    flatmass_cuts.emplace_back(get_flat_mass(content.at(pt_range)));
  }
  return {pt_bins, flatmass_cuts};
}

GN2XHandler::FlatMass GN2XHandler::get_flat_mass(nlohmann::json j)
{
  std::vector<int> mass;
  auto mass_array = j.at("mass");
  for (const auto& element : mass_array) {
    // Skip non-numeric elements like "inf"
    if (element.is_number_integer()) {
      mass.push_back(element.get<int>());
    }
  }
  const std::vector<float> cutvalues = j.at("cutvalues").get<std::vector<float>>();
  return FlatMass {mass, cutvalues};
}

float GN2XHandler::get_cut(const MassPtCuts& cuts, float mass, float pt)
{
  const std::vector<int>& pt_bins = cuts.first;
  int pt_bin_index = find_pt_bin(pt_bins, pt);

  const auto& flatmass = cuts.second[pt_bin_index];
  const std::vector<int>& mass_bins = flatmass.mass;
  int mass_bin_index = find_mass_bin(mass_bins, mass);

  return flatmass.cutvalues[mass_bin_index];
}

int GN2XHandler::find_pt_bin(const std::vector<int>& bin_edges, float pt)
{
  // Assumes last pt bin is [bin_edges.back(), infinity)

  auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), pt);
  int index = std::distance(bin_edges.begin(), it) - 1;

  // Clamp index to the valid range [0, number_of_bins - 1]
  if (index < 0) {
    spdlog::warn("p_T below valid range: {} < {}", pt, bin_edges[0]);
    index = 0;
  } else if (index >= bin_edges.size() - 1) {
    // No warning, it is ok and expected
    index = bin_edges.size() - 1;
  }

  return index;
}

int GN2XHandler::find_mass_bin(const std::vector<int>& bin_edges, float mass)
{
  int start = bin_edges[0];
  int end = bin_edges[1];
  // Assuming uniform bins, use the second edge to find step
  int step = end - start;

  int num_bins = static_cast<int>(bin_edges.size() - 1);
  int index = static_cast<int>(std::floor((mass - start) / step));

  // Clamp index to the valid range [0, num_bins - 1]
  if (index < 0) {
    // this is fine
    index = 0;
  } else if (index >= num_bins) {
    spdlog::warn("mass below valid range: {} < {}", mass, bin_edges.back());
    index = num_bins - 1;
  }

  return index;
}

ROOT::RDF::RNode GN2XHandler::define_pass(ROOT::RDF::RNode df, const std::string& dhbb, const std::string& ljet_m, const std::string& ljet_pt)
{
  const float scale_factor = 1;
  for(const auto& [wp, cuts] : all_cuts_) {
    std::string new_col_name = "pass_" + wp;

    auto apply_cut = [this, &cuts, scale_factor](float dhbb, float ljet_m, float ljet_pt) {
      float scaled_m = ljet_m * scale_factor;
      float scaled_pt = ljet_pt * scale_factor;
      float cut_value = this->get_cut(cuts, scaled_m, scaled_pt);
      bool pass_dhbb = dhbb > cut_value;
      return pass_dhbb;
    };
    df = df.Define(new_col_name, apply_cut, {dhbb, ljet_m, ljet_pt});
  }
	return df;
}

ROOT::RDF::RNode GN2XHandler::define_pass_rvec(ROOT::RDF::RNode df, const std::string& dhbb, const std::string& ljet_m, const std::string& ljet_pt)
{
	const float scale_factor = 1.0;
	for (const auto& [wp, cuts] : all_cuts_) {
		std::string new_col_name = "pass_" + wp;

		auto apply_cut_vectorized = [this, &cuts, scale_factor](
				const ROOT::RVec<float>& dhbb_vec,
				const ROOT::RVec<float>& ljet_m_vec,
				const ROOT::RVec<float>& ljet_pt_vec)
		{
			size_t n_elements = dhbb_vec.size();
			ROOT::RVec<bool> pass_flags(n_elements);

			for (size_t i = 0; i < n_elements; ++i) {
				float scaled_m = ljet_m_vec[i] * scale_factor;
				float scaled_pt = ljet_pt_vec[i] * scale_factor;
				float cut_value = this->get_cut(cuts, scaled_m, scaled_pt);
				pass_flags[i] = (dhbb_vec[i] > cut_value);
			}
			return pass_flags;
		};
		df = df.Define(new_col_name, apply_cut_vectorized, {dhbb, ljet_m, ljet_pt});
	}
	return df;
}
