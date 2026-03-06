set(BJES_SOURCES
    main_bjes.cpp
    config_bjes.cpp
    config_bjes.h
    bjes.cpp
    bjes.h
)

add_executable(run_bjes ${BJES_SOURCES})
target_include_directories(run_bjes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(run_bjes PRIVATE runanalysis_core)

set_target_properties(run_bjes PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

