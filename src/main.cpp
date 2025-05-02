#include <string>
#include <cxxopts.hpp>
#include "xbby_calib.h"

int main(int argc, char* argv[])
{
  cxxopts::Options options("run_analysis", "Create histograms for Zbb+gamma calibration");

  options.add_options()
    ("c,config", "Path to JSON config file", cxxopts::value<std::string>())
    ("i,input", "Input folder name", cxxopts::value<std::string>())
    ("o,output", "Output folder name", cxxopts::value<std::string>()->default_value("output"))
    ("m,multi-threading", "Enable multi-threading", cxxopts::value<bool>()->default_value("true"))
    ("h,help", "Print usage");

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  if (!(result.count("input") && result.count("config"))) {
    std::cerr << "Error: config and input arguments are required!\n";
    std::cerr << options.help() << std::endl;
    return 1;
  }

  const auto config_file = result["config"].as<std::string>();
  const auto input_folder = result["input"].as<std::string>();
  const auto multi_threading = result["multi-threading"].as<bool>();
  auto output_folder = result["output"].as<std::string>();
  auto analysis = xbbycalib::Analysis(config_file, input_folder, output_folder, multi_threading);
  analysis.run();

  return 0;
}
