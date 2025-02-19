#pragma once

#include "args_def.hpp"
#include "cte/string.hpp"
#include "cxxopts.hpp"
#include <string>

#include "args_generated.hpp"

struct Args : ArgsGenerated {
  Args() noexcept;
  Args(int argc, char *argv[]) noexcept;

private:
  void parse_options(int argc, char *argv[]) noexcept;
};
