include_directories(${PROJECT_SOURCE_DIR}/src)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors -Wconversion -Wsign-conversion")
if (USE_LIBCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

option(
    TEST_WITH_TSAN
    "By default, address and undefined sanitizers are used for the tests. If
    this is ON, the thread sanitizer is used instead"
    OFF
)

function(create_tests suffix)
  add_library(
      draco_test_${ARGV0} STATIC EXCLUDE_FROM_ALL
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/draco.cpp
      ${PROJECT_SOURCE_DIR}/src/graph_cut.cpp
      ${PROJECT_SOURCE_DIR}/src/kolmogorov_smirnov.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_writer.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/ptba.cpp
      ${PROJECT_SOURCE_DIR}/src/ptba.cpp
      ${PROJECT_SOURCE_DIR}/src/results/analysis.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/windows_merger.cpp
      ${PROJECT_SOURCE_DIR}/src/clusters_replicates.cpp
      ${PROJECT_SOURCE_DIR}/src/compact_ringmap.cpp
      ${PROJECT_SOURCE_DIR}/src/expectation_maximization.cpp
      ${PROJECT_SOURCE_DIR}/src/compact_ringmap.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/reassignment.cpp
      ${PROJECT_SOURCE_DIR}/src/weighted_clusters.cpp
  )
  target_compile_options(draco_test_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(draco_test_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(draco_test_${ARGV0} ${ARGN} ${ARMADILLO_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} fmt::fmt)
  if (NOT dlib_USES_PKGCONFIG)
    target_link_libraries(
      draco_test_${ARGV0}
      dlib::dlib
    )
  else()
    target_link_libraries(
      draco_test_${ARGV0}
      ${dlib_LIBRARIES}
    )
  endif()

  set(test_names
    ringmap_base
    ringmap_shuffle
    ringmap_concat
    ringmap_data
    graph_cut
    weighted_clusters
    matching_indices
    blocking_queue
    windows_merger
    windows_merger_windows
    windows_merger_cache_indices
    ringmap_window
    weibull_fitter
    mutation_map_writer
    draco
    kmeans
    clusters_replicates
    compact_ringmap
    expectation_maximization
    reassignment
  )
  set(ringmap_base_args ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  set(ringmap_shuffle_args ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  set(ringmap_concat_args ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  set(ringmap_data_args ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  set(windows_merger_args ${CMAKE_CURRENT_SOURCE_DIR}/windows_merger_serialized_data.txt)
  set(mutation_map_writer_args ${CMAKE_CURRENT_SOURCE_DIR}/test.mm)

  foreach(test_name IN LISTS test_names)
    set(test_full_name test_${test_name}_${ARGV0})
    add_executable(${test_full_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
    target_compile_options(${test_full_name} PRIVATE ${ARGN})
    target_include_directories(${test_full_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
    target_link_libraries(${test_full_name} ${ARGN} draco_test_${ARGV0})
    add_test(${test_full_name} ${test_full_name} ${${test_name}_args})
    set_tests_properties(${test_full_name} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
    add_dependencies(${test_full_name} run_args_generate)
    add_dependencies(check ${test_full_name})
  endforeach()
endfunction()

if (TEST_WITH_TSAN)
    create_tests(tsan -fsanitize=thread)
else()
    create_tests(aubsan -fsanitize=address,undefined)
endif()
