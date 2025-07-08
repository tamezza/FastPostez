#include <string>
#include <cxxopts.hpp>

#include "bjes.h"
#include "config_bjes.h"

int main(int argc, char* argv[])
{
  cxxopts::Options options("run_analysis", "Create histograms and ntuples for Zbb+gamma calibration");

  options.add_options()
    ("c,config", "Path to JSON config file", cxxopts::value<std::string>())
    ("i,input", "Input folder name", cxxopts::value<std::string>())
    ("o,output", "Output folder name", cxxopts::value<std::string>()->default_value("output"))
    ("n,nthreads", "Number of threads", cxxopts::value<int>()->default_value("1"))
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
  const auto n_threads = result["nthreads"].as<int>();
  auto output_folder = result["output"].as<std::string>();

  const auto config = config_bjes::load_config(config_file);
  auto analysis = bjes::Analysis(config, input_folder, output_folder, n_threads);
  analysis.run();

  return 0;
}
