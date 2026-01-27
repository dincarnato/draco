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

Analysis::Analysis(Args const &args)
    : args(&args), jsonStream(fs::path(args.output_filename())),
      streamer([this] { streamerLoop(); }) {
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

void Analysis::initStream() {
  if (std::exchange(started, true))
    return;

  jsonStream << '{';
  if (args != nullptr) {
    jsonify(jsonStream, "params", *args) << ',';
  }
  jsonify(jsonStream, "transcripts") << ":[";
}

void Analysis::addTranscript(Transcript &&transcript) {
  assert(not stop.load(std::memory_order_acquire));

  std::lock_guard lock(queueMutex);
  transcripts_.push(std::move(transcript));
  queueCv.notify_one();
}

std::queue<Transcript> const &Analysis::transcripts() const noexcept {
  return transcripts_;
}

void Analysis::streamerLoop() noexcept {
  while (not stop.load(std::memory_order_acquire)) {
    auto transcript = [&] {
      std::unique_lock lock(queueMutex);
      if (transcripts_.empty())
        queueCv.wait(lock);

      if (transcripts_.empty())
        return std::optional<Transcript>{};

      std::optional<Transcript> transcript = std::move(transcripts_.front());
      transcripts_.pop();

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
