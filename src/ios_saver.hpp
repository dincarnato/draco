#pragma once

#include <ios>

struct IosSaver {
  IosSaver(std::ios_base& stream) noexcept
      : stream(&stream), flags(stream.flags()) {}

  ~IosSaver() noexcept { stream->setf(flags); }

private:
  std::ios_base* stream;
  std::ios_base::fmtflags flags;
};
