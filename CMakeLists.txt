cmake_minimum_required(VERSION 3.10)
project(RunAnalysis CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core RIO Net ROOTDataFrame Hist Tree)
include(${ROOT_USE_FILE})

include(FetchContent)
FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.2.1
)
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.15.2
)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(cxxopts spdlog json)

set(COMMON_SOURCES
    src/xbby_calib.cpp
    src/xbby_calib.h
    src/input_handler.cpp
    src/input_handler.h
    src/weights.cpp
    src/weights.h
    src/output_manager.cpp
    src/output_manager.h
    src/config.cpp
    src/config.h
    src/config_bjes.cpp
    src/config_bjes.h
    src/gn2x.cpp
    src/gn2x.h
    src/bjes.cpp
    src/bjes.h
)

add_library(runanalysis_common STATIC ${COMMON_SOURCES})

target_link_libraries(runanalysis_common PUBLIC
    cxxopts::cxxopts
    spdlog::spdlog
    nlohmann_json::nlohmann_json
    ROOT::Core
    ROOT::RIO
    ROOT::Net
    ROOT::ROOTDataFrame
    ROOT::Hist
    ROOT::Tree
)

target_include_directories(runanalysis_common PUBLIC src)

add_executable(run_xbby src/main_xbby_calib.cpp)
add_executable(run_bjes src/main_bjes.cpp)

target_link_libraries(run_xbby PRIVATE runanalysis_common)
target_link_libraries(run_bjes PRIVATE runanalysis_common)

set_target_properties(run_xbby run_bjes PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
