#include "args.hpp"

#include <iostream>

Args::Args(int argc, char *argv[]) noexcept { parse_options(argc, argv); }

template <typename Arg>
void add_arg_to_opts(cxxopts::OptionAdder &opt_adder, Arg const &arg) {
  auto description = [&] {
    if constexpr (Arg::default_value_type::is_available) {
      if constexpr (not std::is_same_v<
                        typename Arg::default_value_type::value_type, bool>) {
        if (arg.get_description().max_size == 0) {
          return arg.get_description()
              .append("(Default: ")
              .append(arg.get_default_string())
              .append(")\0");
        } else {
          return arg.get_description()
              .append(" (Default: ")
              .append(arg.get_default_string())
              .append(")");
        }
      } else {
        return arg.get_description();
      }
    } else {
      return arg.get_description();
    }
  }();
  opt_adder(arg.get_parameter_name().c_str(), description.c_str(),
            cxxopts::value<typename Arg::arg_type>());
}

template <typename Group, std::size_t... Idx>
void add_args_to_opts(cxxopts::Options &opts, Group const &group,
                      std::index_sequence<Idx...>) {
  auto opt_adder = opts.add_options(group.description.c_str());
  (add_arg_to_opts(opt_adder, std::get<Idx>(group.args)), ...);

  if (group.description.is_empty()) {
    opt_adder("help", "Displays this help");
  }
}

template <typename Group>
void add_group_to_opts(cxxopts::Options &opts, Group const &group) {
  add_args_to_opts(
      opts, group,
      std::make_index_sequence<std::tuple_size_v<typename Group::args_type>>());
}

template <std::size_t... Idx>
void add_groups_to_opts(cxxopts::Options &opts, std::index_sequence<Idx...>) {
  (add_group_to_opts(opts, std::get<Idx>(args::opts.groups)), ...);
}

void add_groups_to_opts(cxxopts::Options &opts) {
  add_groups_to_opts(
      opts,
      std::make_index_sequence<
          std::tuple_size_v<typename decltype(args::opts)::groups_type>>());
}

template <std::size_t... Idx>
void print_help(cxxopts::Options const &opts, std::index_sequence<Idx...>) {
  std::cout << opts.help(
      {std::get<Idx>(args::opts.groups).description.c_str()...});
}

void print_help(cxxopts::Options const &opts) {
  print_help(
      opts,
      std::make_index_sequence<
          std::tuple_size_v<typename decltype(args::opts)::groups_type>>());
}

template <typename Arg>
void check_arguments_for_arg(cxxopts::ParseResult const &results,
                             Arg const &arg) {
  if constexpr (not Arg::default_value_type::is_available) {
    if (arg.is_mandatory()) {
      auto const &parameter_name = arg.get_parameter_name();
      if (not results.count(parameter_name.c_str())) {
        std::cout << "[!] Error: no argument provided for parameter '"
                  << parameter_name.c_str() << "'\n";
        std::exit(2);
      }
    }
  }
}

template <typename Group, std::size_t... Idx>
void check_arguments_for_group(cxxopts::ParseResult const &results,
                               Group const &group,
                               std::index_sequence<Idx...>) {
  (check_arguments_for_arg(results, std::get<Idx>(group.args)), ...);
}

template <typename Group>
void check_arguments_for_group(cxxopts::ParseResult const &results,
                               Group const &group) {
  check_arguments_for_group(
      results, group,
      std::make_index_sequence<std::tuple_size_v<typename Group::args_type>>());
}

template <std::size_t... Idx>
void check_arguments(cxxopts::ParseResult const &results,
                     std::index_sequence<Idx...>) {
  (check_arguments_for_group(results, std::get<Idx>(args::opts.groups)), ...);
}

void check_arguments(cxxopts::ParseResult const &results) {
  check_arguments(
      results,
      std::make_index_sequence<
          std::tuple_size_v<typename decltype(args::opts)::groups_type>>());
}

void Args::parse_options(int argc, char *argv[]) noexcept {
  cxxopts::Options arg_opts(argv[0], args::opts.description.c_str());

  add_groups_to_opts(arg_opts);

  auto const args_result = [&] {
    try {
      return arg_opts.parse(argc, argv);
    } catch (std::exception &e) {
      std::cerr << "[!] Error: " << e.what() << '\n';
      std::exit(2);
    }
  }();

  if (args_result.count("help")) {
    print_help(arg_opts);
    std::exit(0);
  }

  check_arguments(args_result);
  this->set_parameters_from_args(args_result);
}
