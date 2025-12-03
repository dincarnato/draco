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
  add_executable(ringmap_base_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_base.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp)
  target_compile_options(ringmap_base_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_base_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_base_${ARGV0} ${ARGN} fmt::fmt ${TBB_LIBRARIES})
  
  add_executable(ringmap_shuffle_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_shuffle.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp)
  target_compile_options(ringmap_shuffle_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_shuffle_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_shuffle_${ARGV0} ${ARGN} fmt::fmt ${TBB_LIBRARIES})
  
  add_executable(ringmap_concat_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_concat.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp)
  target_compile_options(ringmap_concat_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_concat_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_concat_${ARGV0} ${ARGN} fmt::fmt ${TBB_LIBRARIES})

  add_executable(ringmap_data_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp)
  target_compile_options(ringmap_data_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_data_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_data_${ARGV0} ${ARGN} fmt::fmt ${TBB_LIBRARIES})
  
  add_executable(graph_cut_${ARGV0} EXCLUDE_FROM_ALL
      graph_cut.cpp
      ${PROJECT_SOURCE_DIR}/src/graph_cut.cpp
      ${PROJECT_SOURCE_DIR}/src/clusters_replicates.cpp
  )
  target_compile_options(graph_cut_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(graph_cut_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(graph_cut_${ARGV0} ${ARGN} ${TBB_LIBRARIES})
  
  add_executable(matching_indices_${ARGV0} EXCLUDE_FROM_ALL
      matching_indices.cpp)
  target_compile_options(matching_indices_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(matching_indices_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(matching_indices_${ARGV0} ${ARGN})
  
  add_executable(weighted_clusters_${ARGV0} EXCLUDE_FROM_ALL
      weighted_clusters.cpp)
  target_compile_options(weighted_clusters_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(weighted_clusters_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(weighted_clusters_${ARGV0} ${ARGN})
  
  add_executable(blocking_queue_${ARGV0} EXCLUDE_FROM_ALL
      blocking_queue.cpp)
  target_compile_options(blocking_queue_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(blocking_queue_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(blocking_queue_${ARGV0} ${ARGN})

  add_executable(windows_merger_${ARGV0} EXCLUDE_FROM_ALL
      windows_merger.cpp
      ${PROJECT_SOURCE_DIR}/src/windows_merger.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp)
  target_compile_options(windows_merger_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(windows_merger_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(windows_merger_${ARGV0} ${ARGN})

  add_executable(windows_merger_windows_${ARGV0} EXCLUDE_FROM_ALL
      windows_merger_windows.cpp)
  target_compile_options(windows_merger_windows_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(windows_merger_windows_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(windows_merger_windows_${ARGV0} ${ARGN})

  add_executable(windows_merger_cache_indices_${ARGV0} EXCLUDE_FROM_ALL
      windows_merger_cache_indices.cpp)
  target_compile_options(windows_merger_cache_indices_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(windows_merger_cache_indices_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(windows_merger_cache_indices_${ARGV0} ${ARGN})

  add_executable(ringmap_window_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_window.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
  )
  target_compile_options(ringmap_window_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_window_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_window_${ARGV0} ${ARGN} fmt::fmt ${TBB_LIBRARIES})

  add_executable(weibull_fitter_${ARGV0} EXCLUDE_FROM_ALL
      weibull_fitter.cpp)
  target_compile_options(weibull_fitter_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(weibull_fitter_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(weibull_fitter_${ARGV0} ${ARGN})

  add_executable(mutation_map_writer_${ARGV0} EXCLUDE_FROM_ALL
      mutation_map_writer.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_writer.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
  )
  target_compile_options(mutation_map_writer_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(mutation_map_writer_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(mutation_map_writer_${ARGV0} ${ARGN})

  add_executable(draco_${ARGV0} EXCLUDE_FROM_ALL
      draco.cpp
      ${PROJECT_SOURCE_DIR}/src/draco.cpp
      ${PROJECT_SOURCE_DIR}/src/args.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/graph_cut.cpp
      ${PROJECT_SOURCE_DIR}/src/windows_merger.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_writer.cpp
      ${PROJECT_SOURCE_DIR}/src/ptba.cpp
      ${PROJECT_SOURCE_DIR}/src/kolmogorov_smirnov.cpp
      ${PROJECT_SOURCE_DIR}/src/results/analysis.cpp
      ${PROJECT_SOURCE_DIR}/src/logger.cpp
      ${PROJECT_SOURCE_DIR}/src/clusters_replicates.cpp
  )
  target_compile_options(draco_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(draco_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(draco_${ARGV0} ${ARGN})

  add_executable(kmeans_${ARGV0} EXCLUDE_FROM_ALL
      kmeans.cpp
  )
  target_compile_options(kmeans_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(kmeans_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(kmeans_${ARGV0} ${ARGN})

  add_executable(clusters_replicates_${ARGV0} EXCLUDE_FROM_ALL
      clusters_replicates.cpp
      ${PROJECT_SOURCE_DIR}/src/clusters_replicates.cpp
  )
  target_compile_options(clusters_replicates_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(clusters_replicates_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(clusters_replicates_${ARGV0} ${ARGN})

  add_test(ringmap_base_${ARGV0} ringmap_base_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(ringmap_shuffle_${ARGV0} ringmap_shuffle_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(ringmap_concat_${ARGV0} ringmap_concat_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(ringmap_data_${ARGV0} ringmap_data_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(graph_cut_${ARGV0} graph_cut_${ARGV0})
  add_test(matching_indices_${ARGV0} matching_indices_${ARGV0})
  add_test(weighted_clusters_${ARGV0} weighted_clusters_${ARGV0})
  add_test(blocking_queue_${ARGV0} blocking_queue_${ARGV0})
  add_test(windows_merger_${ARGV0} windows_merger_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/windows_merger_serialized_data.txt)
  add_test(windows_merger_windows_${ARGV0} windows_merger_windows_${ARGV0})
  add_test(windows_merger_cache_indices_${ARGV0} windows_merger_cache_indices_${ARGV0})
  add_test(ringmap_window_${ARGV0} ringmap_window_${ARGV0})
  add_test(weibull_fitter_${ARGV0} weibull_fitter_${ARGV0})
  add_test(mutation_map_writer_${ARGV0} mutation_map_writer_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/test.mm)
  add_test(draco_${ARGV0} draco_${ARGV0})
  add_test(kmeans_${ARGV0} kmeans_${ARGV0})
  add_test(clusters_replicates_${ARGV0} clusters_replicates_${ARGV0})

  set_tests_properties(ringmap_base_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_shuffle_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_concat_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_data_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(graph_cut_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(matching_indices_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(weighted_clusters_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(blocking_queue_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_windows_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_cache_indices_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_window_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(weibull_fitter_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(mutation_map_writer_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(draco_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(kmeans_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(clusters_replicates_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  
  target_link_libraries(windows_merger_${ARGV0} ${TBB_LIBRARIES} fmt::fmt)
  target_link_libraries(ringmap_base_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_data_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(graph_cut_${ARGV0} ${ARMADILLO_LIBRARIES} ${TBB_LIBRARIES} fmt::fmt)
  target_link_libraries(ringmap_shuffle_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_concat_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_window_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(weibull_fitter_${ARGV0} dlib::dlib)
  target_link_libraries(draco_${ARGV0} ${ARMADILLO_LIBRARIES} ${TBB_LIBRARIES} ${Boost_LIBRARIES} fmt::fmt)
  if (NOT dlib_USES_PKGCONFIG)
    target_link_libraries(
      draco_${ARGV0}
      dlib::dlib
    )
  else()
    target_link_libraries(
      draco_${ARGV0}
      ${dlib_LIBRARIES}
    )
  endif()
  target_link_libraries(kmeans_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(clusters_replicates_${ARGV0} ${ARMADILLO_LIBRARIES} fmt::fmt)


  add_dependencies(ringmap_base_${ARGV0} run_args_generate)
  add_dependencies(ringmap_shuffle_${ARGV0} run_args_generate)
  add_dependencies(ringmap_concat_${ARGV0} run_args_generate)
  add_dependencies(ringmap_data_${ARGV0} run_args_generate)
  add_dependencies(graph_cut_${ARGV0} run_args_generate)
  add_dependencies(matching_indices_${ARGV0} run_args_generate)
  add_dependencies(weighted_clusters_${ARGV0} run_args_generate)
  add_dependencies(blocking_queue_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_windows_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_cache_indices_${ARGV0} run_args_generate)
  add_dependencies(ringmap_window_${ARGV0} run_args_generate)
  add_dependencies(weibull_fitter_${ARGV0} run_args_generate)
  add_dependencies(mutation_map_writer_${ARGV0} run_args_generate)
  add_dependencies(draco_${ARGV0} run_args_generate)
  add_dependencies(clusters_replicates_${ARGV0} run_args_generate)
  
  add_dependencies(check
    ringmap_base_${ARGV0}
    ringmap_shuffle_${ARGV0}
    ringmap_concat_${ARGV0}
    ringmap_data_${ARGV0}
    graph_cut_${ARGV0}
    matching_indices_${ARGV0}
    weighted_clusters_${ARGV0}
    blocking_queue_${ARGV0}
    windows_merger_${ARGV0}
    windows_merger_windows_${ARGV0}
    windows_merger_cache_indices_${ARGV0}
    ringmap_window_${ARGV0}
    weibull_fitter_${ARGV0}
    mutation_map_writer_${ARGV0}
    draco_${ARGV0}
    kmeans_${ARGV0}
    clusters_replicates_${ARGV0}
  )
endfunction()

if (TEST_WITH_TSAN)
    create_tests(tsan -fsanitize=thread)
else()
    create_tests(aubsan -fsanitize=address,undefined)
endif()
