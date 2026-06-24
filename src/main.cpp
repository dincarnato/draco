#include "args.hpp"
#include "draco.hpp"
#include "logger.hpp"
#include "mutation_map.hpp"
#include "parallel/blocking_queue.hpp"
#include "results/analysis.hpp"
#include "ringmap_data.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <optional>
#include <ranges>
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

  {
    auto window_sizes = args.window_size();
    auto window_sizes_size = std::size(window_sizes);
    auto transcript_fractions = args.window_size_fraction_transcript_size();
    auto transcript_fractions_size = std::size(transcript_fractions);
    auto using_fractions =
        transcript_fractions_size > 1 or transcript_fractions[0] != 0;
    if (using_fractions and window_sizes_size > 1) {
      bail("Invalid parameters: cannot specify both --winLen and "
           "--winLenFracRnaLen");
    }

    if (using_fractions) {
      if (std::ranges::any_of(transcript_fractions, [](auto fraction) {
            return fraction <= 0. or fraction > 1.;
          })) {
        bail("Invalid parameters: --minWindowBases can only be used with "
             "values between 0 (excluded) and 1 (included)");
      }
    } else {
      for (auto window_size : window_sizes) {
        if (window_size > 1 and window_size < args.min_bases_size()) {
          bail("Invalid parameters: --winLen ({}) < --minWindowBases ({})",
               window_size, args.min_bases_size());
        }
      }
    }

    auto usable_windows_sizes =
        using_fractions ? transcript_fractions_size : window_sizes_size;
    auto window_shifts = args.window_shift();
    auto window_size_fraction_transcript_sizes =
        args.window_size_fraction_transcript_size();
    auto window_shifts_size = std::size(window_shifts);
    if (window_shifts_size > 1 and window_shifts_size != usable_windows_sizes) {
      bail("Invalid parameters: multiple values given for winOffset ({}), but "
           "this is not coherent with the number of window lengths ({})",
           window_shifts_size, usable_windows_sizes);
    }
  }

  if (args.max_clusters() < 1) {
    bail("--maxClusters must be at least 1");
  }

  if (not args.assignments_dump_directory().empty()) {
    logger::debug("Creating directory {}", args.assignments_dump_directory());
    std::filesystem::create_directories(args.assignments_dump_directory());
  }

  logger::info("Starting DRACO analysis. This might take a while.");

  // Disabling OMP, we parallelize a higher level
  omp_set_num_threads(1);

  std::vector<MutationMap> mutation_maps;
  for (auto const &mm_filename : args.mm_filenames()) {
    mutation_maps.emplace_back(mm_filename);
  }
  auto analysisResult = ([&] {
    try {
      return results::Analysis(args);
    } catch (std::exception &err) {
      bail("unable to create an analysis result: {}", err.what());
    }
  })();

  analysisResult.filenames = args.mm_filenames();

  if (args.create_eigengaps_plots()) {
    fs::create_directory(fs::path(args.eigengaps_plots_root_dir()));
  }

  std::vector<parallel::blocking_queue<std::tuple<
      MutationMapTranscript, RingmapData, std::optional<std::uint16_t>>>>
      queues;
  std::vector<std::thread> readers;
  queues.reserve(std::size(mutation_maps));
  readers.reserve(std::size(mutation_maps));
  for (std::size_t index = 0; index < std::size(mutation_maps); ++index) {
    queues.emplace_back(10);
    readers.emplace_back([&, index] {
      RingmapData::enqueueRingmapsFromMutationMap(mutation_maps[index],
                                                  queues[index], args);
    });
  }

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
      n_processors =
          static_cast<unsigned>(std::max(tbb::info::default_concurrency(), 1));
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

  using popped_data_t =
      std::vector<std::optional<std::tuple<MutationMapTranscript, RingmapData,
                                           std::optional<std::uint16_t>>>>;
  tbb::parallel_pipeline(
      max_allowed_parallelism,
      tbb::make_filter<void, popped_data_t>(
          tbb::filter_mode::serial_in_order,
          [&](tbb::flow_control &flow_control) {
            auto popped_data =
                queues |
                std::views::transform([](auto &queue) { return queue.pop(); }) |
                std::views::as_rvalue | std::ranges::to<std::vector>();
            if (std::ranges::any_of(popped_data,
                                    [](auto const &single_popped_data) {
                                      return not single_popped_data;
                                    })) {
              flow_control.stop();
              return std::vector<
                  std::optional<std::tuple<MutationMapTranscript, RingmapData,
                                           std::optional<std::uint16_t>>>>();
            } else {
              return popped_data;
            }
          }),
      tbb::make_filter<popped_data_t, void>(
          tbb::filter_mode::parallel, [&](popped_data_t &&popped_data) {
            std::vector<MutationMapTranscript const *> transcripts;
            std::vector<RingmapData *> ringmaps_data;
            std::optional<std::uint16_t> forced_window_length;

            if (!popped_data.empty()) {
              transcripts.push_back(&std::get<0>(*popped_data[0]));
              ringmaps_data.push_back(&std::get<1>(*popped_data[0]));
              forced_window_length = std::get<2>(*popped_data[0]);
            }

            for (auto &single_popped_data : popped_data | std::views::drop(1)) {
              transcripts.push_back(&std::get<0>(*single_popped_data));
              ringmaps_data.push_back(&std::get<1>(*single_popped_data));

              if (forced_window_length != std::get<2>(*single_popped_data)) {
                assert(not transcripts.empty());
                bail("Inconsitent forced_window_length on transcript {}",
                     transcripts[0]->getId());
              }
            }

            assert(not transcripts.empty());
            assert(not ringmaps_data.empty());
            auto const &first_transcript = *transcripts[0];
            if (std::ranges::any_of(transcripts | std::views::drop(1),
                                    [&](auto const &transcript) {
                                      return transcript->getId() !=
                                                 first_transcript.getId() or
                                             transcript->getSequence() !=
                                                 first_transcript.getSequence();
                                    })) {
              bail("Expected ordered transcripts between files, discrepancies "
                   "found on transcript {}, exiting.",
                   first_transcript.getId());
            }

            handle_transcripts(transcripts, ringmaps_data, forced_window_length,
                               analysisResult, args, raw_n_clusters_stream,
                               raw_n_clusters_stream_mutex);
          }));

  for (auto &reader : readers) {
    reader.join();
  }

  logger::info("All done.");
}
