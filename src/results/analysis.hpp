#pragma once

#include "transcript.hpp"

#include <atomic>
#include <condition_variable>
#include <fstream>
#include <mutex>
#include <queue>
#include <string>
#include <string_view>
#include <thread>

namespace results {

struct Analysis final {
  Analysis() = default;
  Analysis(std::string_view jsonFilename);
  ~Analysis() noexcept;

  void addTranscript(Transcript &&transcript);
  std::queue<Transcript> const &transcripts() const noexcept;

  std::vector<std::string> filenames;

private:
  void streamerLoop() noexcept;
  void initStream();

  bool started = false;
  bool writtenFirstTranscript = false;
  std::atomic_bool stop{false};
  mutable std::mutex queueMutex;
  std::condition_variable queueCv;
  std::queue<Transcript> transcripts_;
  std::ofstream jsonStream;
  std::thread streamer;
};

} // namespace results
