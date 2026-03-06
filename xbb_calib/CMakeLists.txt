set(XBB_SOURCES
    main.cpp
    config.cpp
    config.h
    xbby_calib.cpp
    xbby_calib.h
)

add_executable(run_xbby ${XBB_SOURCES})
target_include_directories(run_xbby PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(run_xbby PRIVATE runanalysis_core)

set_target_properties(run_xbby PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

