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

namespace args {
struct Args;
} // namespace args

namespace results {

struct Analysis final {
  Analysis() = default;
  Analysis(Args const &args);
  ~Analysis() noexcept;

  void addTranscript(Transcript &&transcript);
  std::queue<Transcript> const &transcripts() const noexcept;

  std::vector<std::string> filenames;

private:
  void streamerLoop() noexcept;
  void initStream();

  Args const *args{};

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
