#pragma once

#include <tuple>
#include <type_traits>

#include "cte/string.hpp"
#include "default_value.hpp"

namespace args {

enum class Optionality {
  Mandatory,
  Optional,
};

template <typename Type, typename Default, std::size_t TypenameSize,
          std::size_t VariableSize, std::size_t ParameterSize,
          std::size_t DescriptionSize>
struct Arg {
  static_assert(is_default_value_v<Default>);
  static_assert(!Default::is_available or
                is_convertible_from_default_value_v<Type, Default>);
  template <typename, typename, std::size_t, std::size_t, std::size_t,
            std::size_t>
  friend struct Arg;
  using arg_type = Type;
  using default_value_type = Default;

  constexpr Arg(cte::string<TypenameSize> type_name,
                cte::string<VariableSize> variable_name) noexcept
      : _type_name(type_name), _variable_name(variable_name),
        _parameter_name(variable_name.replace('_', '-')), _description(),
        _optionality(Optionality::Mandatory),
        _default_value(MAKE_NO_DEFAULT_VALUE) {}

  template <std::size_t _ParameterSize>
  constexpr auto
  parameter_name(cte::string<_ParameterSize> parameter_name) noexcept {
    return Arg<Type, Default, TypenameSize, VariableSize, _ParameterSize,
               DescriptionSize>{
        std::move(_type_name),     std::move(_variable_name),
        std::move(parameter_name), std::move(_description),
        std::move(_optionality),   std::move(_default_value)};
  }

  template <std::size_t ArgSize>
  constexpr auto
  parameter_name(char const (&parameter_name)[ArgSize]) noexcept {
    return Arg<Type, Default, TypenameSize, VariableSize, ArgSize - 1,
               DescriptionSize>{
        std::move(_type_name),       std::move(_variable_name),
        cte::string(parameter_name), std::move(_description),
        std::move(_optionality),     std::move(_default_value)};
  }

  template <std::size_t _DescriptionSize>
  constexpr auto
  description(cte::string<_DescriptionSize> description) noexcept {
    return Arg<Type, Default, TypenameSize, VariableSize, ParameterSize,
               _DescriptionSize>{
        std::move(_type_name),      std::move(_variable_name),
        std::move(_parameter_name), std::move(description),
        std::move(_optionality),    std::move(_default_value)};
  }

  template <std::size_t ArgSize>
  constexpr auto description(char const (&description)[ArgSize]) noexcept {
    return Arg<Type, Default, TypenameSize, VariableSize, ParameterSize,
               ArgSize - 1>{
        std::move(_type_name),        std::move(_variable_name),
        cte::string(_parameter_name), std::move(description),
        std::move(_optionality),      std::move(_default_value)};
  }

  constexpr auto optional() noexcept {
    return Arg{std::move(_type_name),      std::move(_variable_name),
               std::move(_parameter_name), std::move(_description),
               Optionality::Optional,      std::move(_default_value)};
  }

  constexpr auto mandatory() noexcept {
    return Arg{std::move(_type_name),      std::move(_variable_name),
               std::move(_parameter_name), std::move(_description),
               Optionality::Mandatory,     std::move(_default_value)};
  }

  constexpr auto optionality(Optionality optionality) noexcept {
    return Arg{std::move(_type_name),      std::move(_variable_name),
               std::move(_parameter_name), std::move(_description),
               std::move(optionality),     std::move(_default_value)};
  }

  template <typename _Default>
  constexpr auto default_value(_Default value) noexcept {
    return Arg<Type, _Default, TypenameSize, VariableSize, ParameterSize,
               DescriptionSize>{
        std::move(_type_name),      std::move(_variable_name),
        std::move(_parameter_name), std::move(_description),
        std::move(_optionality),    value};
  }

  constexpr cte::string<TypenameSize> const &get_type_name() const noexcept {
    return _type_name;
  }

  constexpr cte::string<VariableSize> const &
  get_variable_name() const noexcept {
    return _variable_name;
  }

  constexpr cte::string<ParameterSize> const &
  get_parameter_name() const noexcept {
    return _parameter_name;
  }

  constexpr cte::string<DescriptionSize> const &
  get_description() const noexcept {
    return _description;
  }

  constexpr bool is_optional() const noexcept {
    return _optionality == Optionality::Optional;
  }

  constexpr bool is_mandatory() const noexcept {
    return _optionality == Optionality::Mandatory;
  }

  constexpr auto get_default_value() const noexcept {
    return _default_value.value();
  }

  constexpr auto get_default_string() const noexcept {
    return _default_value.into_string();
  }

protected:
  constexpr Arg(cte::string<TypenameSize> type_name,
                cte::string<VariableSize> variable_name,
                cte::string<ParameterSize> parameter_name,
                cte::string<DescriptionSize> description,
                Optionality optionality, Default default_value) noexcept
      : _type_name(std::move(type_name)),
        _variable_name(std::move(variable_name)),
        _parameter_name(std::move(parameter_name)),
        _description(std::move(description)),
        _optionality(std::move(optionality)),
        _default_value(std::move(default_value)) {}

  cte::string<TypenameSize> _type_name;
  cte::string<VariableSize> _variable_name;
  cte::string<ParameterSize> _parameter_name;
  cte::string<DescriptionSize> _description;
  Optionality _optionality;
  Default _default_value;
};

template <typename T> struct is_arg : std::false_type {};

template <typename Type, typename Default, std::size_t TypenameSize,
          std::size_t VariableSize, std::size_t ParameterSize,
          std::size_t DescriptionSize>
struct is_arg<Arg<Type, Default, TypenameSize, VariableSize, ParameterSize,
                  DescriptionSize>> : std::true_type {};

template <typename T> constexpr bool is_arg_v = is_arg<T>::value;

template <std::size_t DescriptionSize, typename... Args> struct Group {
  static_assert(std::conjunction_v<is_arg<Args>...>);
  static constexpr std::size_t description_size = DescriptionSize;
  using args_type = std::tuple<Args...>;

  template <typename... _Args>
  constexpr Group(cte::string<DescriptionSize> description,
                  _Args &&...args) noexcept
      : description(std::move(description)),
        args(std::forward<_Args>(args)...) {}

  template <typename... _Args>
  constexpr Group(char const (&description)[DescriptionSize + 1],
                  _Args &&...args) noexcept
      : description(cte::string(description)),
        args(std::forward<_Args>(args)...) {}

  cte::string<DescriptionSize> description;
  args_type args;
};

template <std::size_t DescriptionSize, typename... Args>
Group(cte::string<DescriptionSize>, Args &&...)
    -> Group<DescriptionSize, std::decay_t<Args>...>;

template <std::size_t StringSize, typename... Args>
Group(char const (&)[StringSize], Args &&...)
    -> Group<StringSize - 1, std::decay_t<Args>...>;

template <typename T> struct is_group : std::false_type {};

template <std::size_t DescriptionSize, typename... Args>
struct is_group<Group<DescriptionSize, Args...>> : std::true_type {};

template <typename T> constexpr bool is_group_v = is_group<T>::value;

template <std::size_t DescriptionSize, typename... Groups> struct Opts {
  static_assert(std::conjunction_v<is_group<Groups>...>);
  static constexpr std::size_t description_size = DescriptionSize;
  using groups_type = std::tuple<Groups...>;

  template <typename... _Groups>
  constexpr Opts(cte::string<DescriptionSize> description,
                 _Groups &&...groups) noexcept
      : description(std::move(description)),
        groups(std::forward<_Groups>(groups)...) {}

  template <typename... _Groups>
  constexpr Opts(char const (&description)[DescriptionSize + 1],
                 _Groups &&...groups) noexcept
      : description(description), groups(std::forward<_Groups>(groups)...) {}

  cte::string<DescriptionSize> description;
  groups_type groups;
};

template <std::size_t DescriptionSize, typename... Groups>
Opts(cte::string<DescriptionSize> description, Groups &&...groups)
    -> Opts<DescriptionSize, Groups...>;

template <std::size_t DescriptionSize, typename... Groups>
Opts(char const (&description)[DescriptionSize], Groups &&...groups)
    -> Opts<DescriptionSize - 1, Groups...>;

} // namespace args

#define ARG(type, name)                                                        \
  ::args::Arg<type, ::DefaultValue<0, ::DefaultValueType::None>,               \
              ::cte::string(#type).size(), ::cte::string(#name).size(),        \
              ::cte::string(#name).size(), ::std::size_t(0)>(                  \
      cte::string(#type), cte::string(#name))

#define DEFAULT_VALUE(x) default_value(MAKE_DEFAULT_VALUE(x))
