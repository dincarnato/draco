#include "args.hpp"
#include "draco.hpp"
#include "logger.hpp"
#include "mutation_map.hpp"
#include "parallel/blocking_queue.hpp"
#include "results/analysis.hpp"
#include "ringmap_data.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <optional>
#include <thread>

#include <armadillo>
#include <tbb/global_control.h>

#include <omp.h>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  auto const args = Args(argc, argv);

  auto &&log_level = args.log_level();
  if (not log_level.empty()) {
    auto level = logger::parse_level(log_level);
    if (level.has_value()) {
      logger::instance.set_level(*level);
    } else {
      logger::warn("Invalid debug level \"{}\"", log_level);
    }
  }

  if (not args.assignments_dump_directory().empty()) {
    logger::debug("Creating directory {}", args.assignments_dump_directory());
    std::filesystem::create_directories(args.assignments_dump_directory());
  }

  std::cout << "\n[+] Starting DRACO analysis. This might take a while...\n";

  // Disabling OMP, we parallelize a higher level
  omp_set_num_threads(1);

  MutationMap mutationMap(args.mm_filename());
  results::Analysis analysisResult(args.output_filename());

  analysisResult.filename = args.mm_filename();

  if (args.create_eigengaps_plots()) {
    fs::create_directory(fs::path(args.eigengaps_plots_root_dir()));
  }

  parallel::blocking_queue<std::pair<MutationMapTranscript, RingmapData>> queue(
      10);
  std::thread reader([&] {
    RingmapData::enqueueRingmapsFromMutationMap(mutationMap, queue, args);
  });
  /*
  std::thread reader([&] {
    auto transcriptIter = std::next(std::begin(mutationMap), 13);
    queue.emplace(*transcriptIter, RingmapData(*transcriptIter, args));
    ++transcriptIter;
    queue.emplace(*transcriptIter, RingmapData(*transcriptIter, args));
    queue.finish();
  });
  */

  const auto max_allowed_parallelism = [&] {
    auto n_processors = args.n_processors();
    if (n_processors == 0) {
      n_processors = tbb::info::default_concurrency();
    }

    return static_cast<std::size_t>(n_processors);
  }();

  tbb::global_control tbb_global_control(
      tbb::global_control::max_allowed_parallelism, max_allowed_parallelism);

  std::optional<std::ofstream> raw_n_clusters_stream;
  std::mutex raw_n_clusters_stream_mutex;
  if (!args.output_raw_n_clusters().empty()) {
    raw_n_clusters_stream =
        std::ofstream(args.output_raw_n_clusters(),
                      std::ios_base::out | std::ios_base::trunc);
  }

  tbb::parallel_pipeline(
      max_allowed_parallelism,
      tbb::make_filter<void, void>(
          tbb::filter_mode::parallel, [&](tbb::flow_control &flow_control) {
            auto poppedData = queue.pop();
            if (not poppedData) {
              flow_control.stop();
              return;
            }

            auto const &transcript = std::get<0>(*poppedData);
            auto &ringmapData = std::get<1>(*poppedData);

            handle_transcript(transcript, ringmapData, analysisResult, args,
                              raw_n_clusters_stream,
                              raw_n_clusters_stream_mutex);
          }));

  reader.join();
  std::cout << "\n[+] All done.\n\n";
}
