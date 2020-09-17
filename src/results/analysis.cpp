#include "analysis.hpp"
#include "transcript.hpp"
#include "window.hpp"

#include "jsonify.hpp"

#include <cassert>
#include <optional>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "Missing filesystem header"
#endif

namespace results {

Analysis::Analysis(std::string_view jsonFilename)
    : jsonStream(fs::path(jsonFilename)), streamer([this] { streamerLoop(); }) {
  if (jsonStream.fail())
    throw std::ofstream::failure("output json file cannot be opened");
}

Analysis::~Analysis() noexcept {
  queueMutex.lock();
  stop = true;
  if (streamer.joinable()) {
    queueCv.notify_one();
    queueMutex.unlock();
    streamer.join();
  }

  if (jsonStream) {
    if (not started)
      initStream();

    jsonStream << "]}";
  }
}

void
Analysis::initStream() {
  if (std::exchange(started, true))
    return;

  jsonStream << '{';
  jsonify(jsonStream, "filename", filename) << ',';
  jsonify(jsonStream, "transcripts") << ":[";
}

void
Analysis::addTranscript(Transcript&& transcript) {
  assert(not stop.load(std::memory_order_relaxed));

  std::lock_guard lock(queueMutex);
  transcripts.push(std::move(transcript));
  queueCv.notify_one();
}

void
Analysis::streamerLoop() noexcept {
  while (not stop.load(std::memory_order_relaxed)) {
    auto transcript = [&] {
      std::unique_lock lock(queueMutex);
      if (transcripts.empty())
        queueCv.wait(lock);

      if (transcripts.empty())
        return std::optional<Transcript>{};

      std::optional<Transcript> transcript = std::move(transcripts.front());
      transcripts.pop();

      return transcript;
    }();

    if (not transcript) {
      if (stop.load(std::memory_order_acquire))
        return;

      continue;
    }

    if (not started)
      initStream();

    if (std::exchange(writtenFirstTranscript, true))
      jsonStream << ',';
    jsonify(jsonStream, *transcript);
  }
}

} // namespace results
