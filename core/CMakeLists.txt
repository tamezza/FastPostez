set(CORE_SOURCES
    input_handler.cpp
    input_handler.h
    output_manager.cpp
    output_manager.h
    weights.cpp
    weights.h
    gn2x.cpp
    gn2x.h
)

add_library(runanalysis_core STATIC ${CORE_SOURCES})
target_include_directories(runanalysis_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(runanalysis_core PUBLIC
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

