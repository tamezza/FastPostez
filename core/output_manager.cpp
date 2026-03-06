#include <string>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <memory>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "output_manager.h"

namespace fs = std::filesystem;

std::string output_manager::create_output_folder(const std::string& output_folder)
{
  const std::string& output_folder_timestamp = output_folder + "_" + output_manager::get_timestamp();
  fs::path base_dir(output_folder_timestamp);
  fs::path histo_dir = base_dir / "histograms";
  fs::path ntuple_dir = base_dir / "ntuples";

  fs::create_directories(histo_dir);
  fs::create_directories(ntuple_dir);

  return output_folder_timestamp;
}

std::string output_manager::get_timestamp()
{
  auto now = std::chrono::system_clock::now();
  auto time_t_now = std::chrono::system_clock::to_time_t(now);
  std::tm tm_now = *std::localtime(&time_t_now);

  std::ostringstream timestamp;
  timestamp << std::put_time(&tm_now, "%Y-%m-%d_%H-%M-%S");
  return timestamp.str();
}

void output_manager::set_logger(const std::string& output_folder)
{
  const std::string log_file = output_folder + "/analysis.log";
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);

  auto logger_ptr = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list{console_sink, file_sink});
  logger_ptr->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");
  spdlog::set_default_logger(logger_ptr);
}
