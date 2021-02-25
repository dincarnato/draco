#include "args_def.hpp"

#include <fstream>
#include <sstream>
#include <string_view>

template <typename Arg>
static void
dump_arg(std::ostream& os, Arg const& arg) {
  auto const type_name = arg.get_type_name().c_str();
  auto const variable_name = arg.get_variable_name().c_str();

  os << "protected:\n"
     << type_name << " _" << variable_name << ";\npublic:\nauto "
     << variable_name << "() const noexcept -> decltype(auto) { return _"
     << variable_name << ";}\n";
}

template <typename Group, std::size_t... Idx>
static void
dump_group_impl(std::ostream& os, Group const& group,
                std::index_sequence<Idx...>) {
  (dump_arg(os, std::get<Idx>(group.args)), ...);
}

template <typename Group>
static void
dump_group(std::ostream& os, Group const& group) {
  dump_group_impl(
      os, group,
      std::make_index_sequence<std::tuple_size_v<typename Group::args_type>>());
}

template <std::size_t... Idx>
static void
dump_opts_impl(std::ostream& os, std::index_sequence<Idx...>) {
  (dump_group(os, std::get<Idx>(args::opts.groups)), ...);
}

static void
dump_opts(std::ostream& os) {
  dump_opts_impl(
      os, std::make_index_sequence<
              std::tuple_size_v<typename decltype(args::opts)::groups_type>>());
}

template <std::size_t GroupIndex, std::size_t ArgIndex, typename Arg>
static void
create_setter_function_for_arg(std::ostream& os, Arg const& arg) {
  auto const variable_name = arg.get_variable_name().c_str();
  auto const parameter_name = arg.get_parameter_name().c_str();
  auto const type_name = arg.get_type_name().c_str();

  os << "if (results.count(\"" << parameter_name << "\")) { _" << variable_name
     << " = results[\"" << parameter_name << "\"].as<" << type_name << ">(); }";

  if constexpr (Arg::default_value_type::is_available) {
    os << " else { _" << variable_name << " = std::get<" << ArgIndex
       << ">(std::get<" << GroupIndex
       << ">(::args::opts.groups).args).get_default_value(); }";
  }
  os << "\n";
}

template <std::size_t GroupIndex, typename Group, std::size_t... Idx>
static void
create_setter_function_for_group(std::ostream& os, Group const& group,
                                 std::index_sequence<Idx...>) {
  (create_setter_function_for_arg<GroupIndex, Idx>(os,
                                                   std::get<Idx>(group.args)),
   ...);
}

template <std::size_t GroupIndex, typename Group>
static void
create_setter_function_for_group(std::ostream& os, Group const& group) {
  create_setter_function_for_group<GroupIndex>(
      os, group,
      std::make_index_sequence<std::tuple_size_v<typename Group::args_type>>());
}

template <std::size_t... Idx>
static void
create_setter_function_for_groups(std::ostream& os,
                                  std::index_sequence<Idx...>) {
  (create_setter_function_for_group<Idx>(os, std::get<Idx>(args::opts.groups)),
   ...);
}

static void
create_setter_function(std::ostream& os) {
  os << "protected:\nvoid set_parameters_from_args(cxxopts::ParseResult const& "
        "results) {\n";
  create_setter_function_for_groups(
      os, std::make_index_sequence<
              std::tuple_size_v<typename decltype(args::opts)::groups_type>>());
  os << "\n}";
}

template <typename Stream> static void generate_on_stream(Stream &stream) {
  stream << "#pragma once\n#include <string>\n\nstruct ArgsGenerated {\n";
  dump_opts(stream);
  create_setter_function(stream);
  stream << "};";
}

int main() {
  constexpr char const *generated_file_name = "args_generated.hpp";
  if (std::ifstream stream(generated_file_name); stream.good()) {
    std::string old_content{std::istreambuf_iterator<char>(stream),
                            std::istreambuf_iterator<char>()};
    std::stringstream new_content_stream;
    generate_on_stream(new_content_stream);

    std::string new_content(new_content_stream.str());
    if (new_content != old_content) {
      std::ofstream out_stream(generated_file_name);
      new_content_stream.seekg(0);
      std::copy(std::istreambuf_iterator<char>(new_content_stream),
                std::istreambuf_iterator<char>(),
                std::ostreambuf_iterator<char>(out_stream));
    }
  } else {
    auto out_stream = std::ofstream(generated_file_name);
    generate_on_stream(out_stream);
  }
}
