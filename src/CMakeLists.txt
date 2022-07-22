set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DARMA_NO_DEBUG")

add_executable(
    draco
    draco.cpp
    ringmap_data.cpp
    ringmap_matrix.cpp
    mutation_map.cpp
    mutation_map_transcript.cpp
    ptba.cpp
    kolmogorov_smirnov.cpp
    partitioner.cpp
    rna_secondary_structure.cpp
    results/analysis.cpp
    graph_cut.cpp
    paired_rna_secondary_structure.cpp
    spectral_partitioner.cpp
    windows_merger.cpp
    args.cpp
    )

add_executable(
    args_generate
    args_generate.cpp
)

target_link_libraries(
    draco
    ${ARMADILLO_LIBRARIES}
    ${TBB_LIBRARIES}
    ${Boost_LIBRARIES}
)

target_include_directories(
    draco
    PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
    PUBLIC
        ${ARMADILLO_INCLUDE_DIRS}
        ${TBB_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
)



if (NOT dlib_USES_PKGCONFIG)
    target_link_libraries(
        draco
        dlib::dlib
    )

    target_include_directories(
        draco
        PUBLIC
            dlib::dlib
    )
else()
    target_link_libraries(
        draco
        ${dlib_LIBRARIES}
    )

    target_include_directories(
        draco
        PUBLIC
            ${dlib_INCLUDE_DIRS}
    )
endif()

add_custom_target(
    run_args_generate
    COMMAND args_generate
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS args_generate
    COMMENT "Generating header file for command line arguments evaluation"
    SOURCES args_generate.cpp
)

add_dependencies(draco run_args_generate)
