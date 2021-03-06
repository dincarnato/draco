include_directories(${PROJECT_SOURCE_DIR}/src)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
set(CMAKE_BUILD_TYPE Debug FORCE)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic-errors -Wconversion -Wsign-conversion -g ${OpenMP_CXX_FLAGS}")
if (USE_LIBCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

function(create_tests suffix)
  add_executable(ringmap_base_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_base.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp)
  target_compile_options(ringmap_base_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_base_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_base_${ARGV0} ${ARGN})
  
  add_executable(ringmap_shuffle_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_shuffle.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp)
  target_compile_options(ringmap_shuffle_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_shuffle_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_shuffle_${ARGV0} ${ARGN})
  
  add_executable(ringmap_concat_${ARGV0} EXCLUDE_FROM_ALL
      ringmap_concat.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_data.cpp
      ${PROJECT_SOURCE_DIR}/src/ringmap_matrix.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map_transcript.cpp
      ${PROJECT_SOURCE_DIR}/src/mutation_map.cpp
      ${PROJECT_SOURCE_DIR}/src/paired_rna_secondary_structure.cpp
      ${PROJECT_SOURCE_DIR}/src/spectral_partitioner.cpp
      ${PROJECT_SOURCE_DIR}/src/rna_secondary_structure.cpp)
  target_compile_options(ringmap_concat_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_concat_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_concat_${ARGV0} ${ARGN})
  
  add_executable(graph_cut_${ARGV0} EXCLUDE_FROM_ALL
      graph_cut.cpp
      ${PROJECT_SOURCE_DIR}/src/graph_cut.cpp)
  target_compile_options(graph_cut_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(graph_cut_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(graph_cut_${ARGV0} ${ARGN})
  
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
      ${PROJECT_SOURCE_DIR}/src/windows_merger.cpp)
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
  )
  target_compile_options(ringmap_window_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(ringmap_window_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(ringmap_window_${ARGV0} ${ARGN})

  add_executable(weibull_fitter_${ARGV0} EXCLUDE_FROM_ALL
      weibull_fitter.cpp)
  target_compile_options(weibull_fitter_${ARGV0} PRIVATE ${ARGN})
  target_include_directories(weibull_fitter_${ARGV0} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src)
  target_link_libraries(weibull_fitter_${ARGV0} ${ARGN})

  add_test(ringmap_base_${ARGV0} ringmap_base_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(ringmap_shuffle_${ARGV0} ringmap_shuffle_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(ringmap_concat_${ARGV0} ringmap_concat_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/ringmap1.txt)
  add_test(graph_cut_${ARGV0} graph_cut_${ARGV0})
  add_test(matching_indices_${ARGV0} matching_indices_${ARGV0})
  add_test(weighted_clusters_${ARGV0} weighted_clusters_${ARGV0})
  add_test(blocking_queue_${ARGV0} blocking_queue_${ARGV0})
  add_test(windows_merger_${ARGV0} windows_merger_${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/windows_merger_serialized_data.txt)
  add_test(windows_merger_windows_${ARGV0} windows_merger_windows_${ARGV0})
  add_test(windows_merger_cache_indices_${ARGV0} windows_merger_cache_indices_${ARGV0})
  add_test(ringmap_window_${ARGV0} ringmap_window_${ARGV0})
  add_test(weibull_fitter_${ARGV0} weibull_fitter_${ARGV0})

  set_tests_properties(ringmap_base_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_shuffle_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_concat_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(graph_cut_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(matching_indices_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(weighted_clusters_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(blocking_queue_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_windows_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(windows_merger_cache_indices_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(ringmap_window_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  set_tests_properties(weibull_fitter_${ARGV0} PROPERTIES ENVIRONMENT TSAN_OPTIONS=halt_on_error=1)
  
  target_link_libraries(windows_merger_${ARGV0} ${TBB_LIBRARIES})
  target_link_libraries(ringmap_base_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(graph_cut_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_shuffle_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_concat_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(ringmap_window_${ARGV0} ${ARMADILLO_LIBRARIES})
  target_link_libraries(weibull_fitter_${ARGV0} ${DLIB_LIBRARIES})

  add_dependencies(ringmap_base_${ARGV0} run_args_generate)
  add_dependencies(ringmap_shuffle_${ARGV0} run_args_generate)
  add_dependencies(ringmap_concat_${ARGV0} run_args_generate)
  add_dependencies(graph_cut_${ARGV0} run_args_generate)
  add_dependencies(matching_indices_${ARGV0} run_args_generate)
  add_dependencies(weighted_clusters_${ARGV0} run_args_generate)
  add_dependencies(blocking_queue_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_windows_${ARGV0} run_args_generate)
  add_dependencies(windows_merger_cache_indices_${ARGV0} run_args_generate)
  add_dependencies(ringmap_window_${ARGV0} run_args_generate)
  add_dependencies(weibull_fitter_${ARGV0} run_args_generate)
  
  add_dependencies(check ringmap_base_${ARGV0} ringmap_shuffle_${ARGV0} ringmap_concat_${ARGV0} graph_cut_${ARGV0} matching_indices_${ARGV0} weighted_clusters_${ARGV0} blocking_queue_${ARGV0} windows_merger_${ARGV0} windows_merger_windows_${ARGV0} windows_merger_cache_indices_${ARGV0} ringmap_window_${ARGV0} weibull_fitter_${ARGV0})
endfunction()

create_tests(aubsan -fsanitize=address,undefined;-O0)
create_tests(tsan -fsanitize=thread;-O0)
create_tests(aubsan_opt -fsanitize=address,undefined;-O3)
create_tests(tsan_opt -fsanitize=thread;-O3)
