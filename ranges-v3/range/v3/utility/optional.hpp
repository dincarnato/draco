/// \file
// Range v3 library
//
//  Copyright Casey Carter 2017
//
//  Use, modification and distribution is subject to the
//  Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// Project home: https://github.com/ericniebler/range-v3
//

#ifndef RANGES_V3_UTILITY_OPTIONAL_HPP
#define RANGES_V3_UTILITY_OPTIONAL_HPP

#include <exception>
#include <initializer_list>
#include <memory>
#include <new>
#include <range/v3/detail/config.hpp>
#include <concepts/concepts.hpp>
#include <range/v3/utility/static_const.hpp>
#include <range/v3/utility/swap.hpp>
#include <range/v3/utility/in_place.hpp>

namespace ranges
{
    template<typename> struct optional;

    struct bad_optional_access
      : std::exception
    {
        virtual const char *what() const noexcept override
        {
            return "bad optional access";
        }
    };

    struct nullopt_t
    {
        struct tag {};
        explicit constexpr nullopt_t(tag) noexcept {}
    };
#if RANGES_CXX_INLINE_VARIABLES >= RANGES_CXX_INLINE_VARIABLES_17
    inline constexpr nullopt_t nullopt{nullopt_t::tag{}};
#else
    /// \cond
    namespace detail
    {
        template<typename>
        struct nullopt_holder
        {
            static constexpr nullopt_t nullopt{nullopt_t::tag{}};
        };
        template<typename T>
        constexpr nullopt_t nullopt_holder<T>::nullopt;
    }
    /// \endcond
    inline namespace
    {
        constexpr auto &nullopt = detail::nullopt_holder<void>::nullopt;
    }
#endif

    /// \cond
    namespace detail
    {
        template<typename = void>
        [[noreturn]] bool throw_bad_optional_access()
        {
            throw bad_optional_access{};
        }

        namespace optional_adl
        {
            template<typename T, bool = std::is_trivially_destructible<T>::value>
            struct optional_storage
            {
                union
                {
                    char dummy_;
                    meta::_t<std::remove_cv<T>> data_;
                };
                bool engaged_;

                constexpr optional_storage() noexcept
                  : optional_storage(tag{},
                        meta::strict_and<
                            detail::is_trivially_default_constructible<T>,
                            detail::is_trivially_copyable<T>>{})
                {}
                CPP_template(typename... Args)(
                    requires Constructible<T, Args...>)
                constexpr explicit optional_storage(in_place_t, Args &&... args)
                    noexcept(std::is_nothrow_constructible<T, Args...>::value)
                  : data_(static_cast<Args &&>(args)...), engaged_{true}
                {}

                constexpr /*c++14*/ void reset() noexcept
                {
                    engaged_ = false;
                }
            private:
                struct tag {};
                constexpr optional_storage(tag, std::false_type) noexcept
                  : dummy_{}, engaged_{false}
                {}
                constexpr optional_storage(tag, std::true_type) noexcept
                  : data_{}, engaged_{false}
                {}
            };

            template<typename T>
            struct optional_storage<T, false>
            {
                union
                {
                    char dummy_;
                    meta::_t<std::remove_cv<T>> data_;
                };
                bool engaged_;

                ~optional_storage()
                {
                    reset();
                }
                constexpr optional_storage() noexcept
                  : dummy_{}, engaged_{false}
                {}
                CPP_template(typename... Args)(
                    requires Constructible<T, Args...>)
                constexpr explicit optional_storage(in_place_t, Args &&... args)
                    noexcept(std::is_nothrow_constructible<T, Args...>::value)
                  : data_(static_cast<Args &&>(args)...), engaged_{true}
                {}
                optional_storage(optional_storage const &) = default;
                optional_storage(optional_storage &&) = default;
                optional_storage &operator=(optional_storage const &) = default;
                optional_storage &operator=(optional_storage &&) = default;

                void reset() noexcept
                {
                    if (engaged_)
                    {
                        data_.~T();
                        engaged_ = false;
                    }
                }
            };

            template<typename T>
            struct optional_base
              : private optional_storage<T>
            {
                using optional_storage<T>::optional_storage;
                using optional_storage<T>::reset;

                constexpr bool has_value() const noexcept
                {
                    return engaged_;
                }
                constexpr /*c++14*/ T &operator*() & noexcept
                {
                    return RANGES_EXPECT(engaged_), data_;
                }
                constexpr T const &operator*() const & noexcept
                {
                    return RANGES_EXPECT(engaged_), data_;
                }
                constexpr /*c++14*/ T &&operator*() && noexcept
                {
                    return RANGES_EXPECT(engaged_), detail::move(data_);
                }
                constexpr /*c++14*/ T const &&operator*() const && noexcept
                {
                    return RANGES_EXPECT(engaged_), detail::move(data_);
                }
                constexpr /*c++14*/ T *operator->() noexcept
                {
                    return RANGES_EXPECT(engaged_), std::addressof(data_);
                }
                constexpr T const *operator->() const noexcept
                {
                    return RANGES_EXPECT(engaged_), std::addressof(data_);
                }
                CPP_member
                constexpr /*c++14*/ auto swap(optional_base &that)
                    noexcept(std::is_nothrow_move_constructible<T>::value &&
                        is_nothrow_swappable<T>::value) ->
                    CPP_ret(void)(
                        requires MoveConstructible<T> && Swappable<T>)
                {
                    constexpr bool can_swap_trivially =
                        !::concepts::adl_swap_detail::is_adl_swappable_<T>::value &&
                        detail::is_trivially_move_constructible<T>::value &&
                        detail::is_trivially_move_assignable<T>::value;

                    swap_(meta::bool_<can_swap_trivially>{}, that);
                }
            protected:
                template<typename... Args>
                auto construct_from(Args &&... args)
                    noexcept(std::is_nothrow_constructible<T, Args...>::value) ->
                    CPP_ret(T &)(
                        requires Constructible<T, Args...>)
                {
                    RANGES_EXPECT(!engaged_);
                    auto const address = static_cast<void *>(std::addressof(data_));
                    ::new (address) T(static_cast<Args &&>(args)...);
                    engaged_ = true;
                    return data_;
                }
                template<typename U>
                constexpr /*c++14*/
                void assign_from(U &&that)
                    noexcept(
                        std::is_nothrow_constructible<T,
                            decltype(*static_cast<U &&>(that))>::value &&
                        std::is_nothrow_assignable<T &,
                            decltype(*static_cast<U &&>(that))>::value)
                {
                    if (!that.has_value())
                        reset();
                    else if (engaged_)
                        data_ = *static_cast<U &&>(that);
                    else
                    {
                        auto const address = static_cast<void *>(std::addressof(data_));
                        ::new (address) T(*static_cast<U &&>(that));
                        engaged_ = true;
                    }
                }
            private:
                constexpr /*c++14*/
                void swap_(std::true_type, optional_base &that) noexcept
                {
                    ranges::swap(
                        static_cast<optional_storage<T> &>(*this),
                        static_cast<optional_storage<T> &>(that));
                }
                constexpr /*c++14*/
                void swap_(std::false_type, optional_base &that)
                    noexcept(std::is_nothrow_move_constructible<T>::value &&
                        is_nothrow_swappable<T>::value)
                {
                    if (that.engaged_ == engaged_)
                    {
                        if (engaged_)
                            ranges::swap(data_, that.data_);
                    }
                    else
                    {
                        auto &src = engaged_ ? *this : that;
                        auto &dst = engaged_ ? that : *this;
                        dst.construct_from(detail::move(src.data_));
                        src.reset();
                    }
                }

                using optional_storage<T>::engaged_;
                using optional_storage<T>::data_;
            };

            template<typename T>
            struct optional_base<T &>
            {
                optional_base() = default;
                template<typename Arg>
                constexpr explicit CPP_ctor(optional_base)(in_place_t, Arg &&arg)(
                    noexcept(true)
                    requires Constructible<T &, Arg>)
                  : ptr_(std::addressof(arg))
                {}
                constexpr bool has_value() const noexcept
                {
                    return ptr_;
                }
                constexpr T &operator*() const noexcept
                {
                    return RANGES_EXPECT(ptr_), *ptr_;
                }
                constexpr T *operator->() const noexcept
                {
                    return RANGES_EXPECT(ptr_), ptr_;
                }
                constexpr /*c++14*/ void reset() noexcept
                {
                    ptr_ = nullptr;
                }
                CPP_member
                constexpr /*c++14*/ auto swap(optional_base &that)
                    noexcept(is_nothrow_swappable<T>::value) ->
                    CPP_ret(void)(
                        requires Swappable<T>)
                {
                    if (ptr_ && that.ptr_)
                        ranges::swap(*ptr_, *that.ptr_);
                    else
                        ranges::swap(ptr_, that.ptr_);
                }
            protected:
                template<typename U>
                constexpr /*c++14*/
                auto construct_from(U &&ref) noexcept ->
                    CPP_ret(T &)(
                        requires ConvertibleTo<U &, T &>)
                {
                    RANGES_EXPECT(!ptr_);
                    ptr_ = std::addressof(ref);
                    return *ptr_;
                }
                template<typename U>
                constexpr /*c++14*/ void assign_from(U &&that)
                {
                    if (ptr_ && that.ptr_)
                        *ptr_ = *that.ptr_;
                    else
                        ptr_ = that.ptr_;
                }
            private:
                T *ptr_ = nullptr;
            };

            template<typename T>
            struct optional_copy
              : optional_base<T>
            {
                optional_copy() = default;
                optional_copy(optional_copy const &that)
                    noexcept(std::is_nothrow_copy_constructible<T>::value)
                {
                    if (that.has_value())
                        this->construct_from(*that);
                }
                optional_copy(optional_copy &&) = default;
                optional_copy &operator=(optional_copy const &) = default;
                optional_copy &operator=(optional_copy &&) = default;

                using optional_base<T>::optional_base;
            };

            template<typename T>
            using copy_construct_layer =
                meta::if_c<
                    std::is_copy_constructible<T>::value &&
                        !detail::is_trivially_copy_constructible<T>::value,
                    optional_copy<T>,
                    optional_base<T>>;

            template<typename T>
            struct optional_move
              : copy_construct_layer<T>
            {
                optional_move() = default;
                optional_move(optional_move const &) = default;
                optional_move(optional_move &&that)
                    noexcept(std::is_nothrow_move_constructible<T>::value)
                {
                    if (that.has_value())
                        this->construct_from(std::move(*that));
                }
                optional_move &operator=(optional_move const &) = default;
                optional_move &operator=(optional_move &&) = default;

                using copy_construct_layer<T>::copy_construct_layer;
            };

            template<typename T>
            using move_construct_layer =
                meta::if_c<
                    std::is_move_constructible<T>::value &&
                        !detail::is_trivially_move_constructible<T>::value,
                    optional_move<T>,
                    copy_construct_layer<T>>;

            template<typename T>
            struct optional_copy_assign
              : move_construct_layer<T>
            {
                optional_copy_assign() = default;
                optional_copy_assign(optional_copy_assign const &) = default;
                optional_copy_assign(optional_copy_assign &&) = default;
                optional_copy_assign &operator=(optional_copy_assign const &that)
                    noexcept(std::is_nothrow_copy_constructible<T>::value &&
                        std::is_nothrow_copy_assignable<T>::value)
                {
                    this->assign_from(that);
                    return *this;
                }
                optional_copy_assign &operator=(optional_copy_assign &&) = default;

                using move_construct_layer<T>::move_construct_layer;
            };

            template<typename T>
            struct deleted_copy_assign
              : move_construct_layer<T>
            {
                deleted_copy_assign() = default;
                deleted_copy_assign(deleted_copy_assign const &) = default;
                deleted_copy_assign(deleted_copy_assign &&) = default;
                deleted_copy_assign &operator=(deleted_copy_assign const &) = delete;
                deleted_copy_assign &operator=(deleted_copy_assign &&) = default;

                using move_construct_layer<T>::move_construct_layer;
            };

            template<typename T>
            using copy_assign_layer = meta::if_c<
                std::is_copy_constructible<T>::value && std::is_copy_assignable<T>::value,
                meta::if_c<
                    std::is_reference<T>::value ||
                      !(detail::is_trivially_copy_constructible<T>::value &&
                        detail::is_trivially_copy_assignable<T>::value),
                    optional_copy_assign<T>,
                    move_construct_layer<T>>,
                deleted_copy_assign<T>>;

            template<typename T>
            struct optional_move_assign
              : copy_assign_layer<T>
            {
                optional_move_assign() = default;
                optional_move_assign(optional_move_assign const &) = default;
                optional_move_assign(optional_move_assign &&) = default;
                optional_move_assign &operator=(optional_move_assign const &) = default;
                optional_move_assign &operator=(optional_move_assign &&that)
                    noexcept(std::is_nothrow_move_constructible<T>::value &&
                        std::is_nothrow_move_assignable<T>::value)
                {
                    this->assign_from(std::move(that));
                    return *this;
                }

                using copy_assign_layer<T>::copy_assign_layer;
            };

            template<typename T>
            struct deleted_move_assign
              : copy_assign_layer<T>
            {
                deleted_move_assign() = default;
                deleted_move_assign(deleted_move_assign const &) = default;
                deleted_move_assign(deleted_move_assign &&) = default;
                deleted_move_assign &operator=(deleted_move_assign const &) = default;
                deleted_move_assign &operator=(deleted_move_assign &&) = delete;

                using copy_assign_layer<T>::copy_assign_layer;
            };

            template<typename T>
            using move_assign_layer = meta::if_c<
                std::is_move_constructible<T>::value && std::is_move_assignable<T>::value,
                meta::if_c<
                    std::is_reference<T>::value ||
                        !(detail::is_trivially_move_constructible<T>::value &&
                          detail::is_trivially_move_assignable<T>::value),
                    optional_move_assign<T>,
                    copy_assign_layer<T>>,
                deleted_move_assign<T>>;
        } // namespace optional_adl
    } // namespace detail
    /// \endcond

    CPP_def
    (
        template(typename U, typename T)
        concept OptionalShouldConvert,
            requires(int)(void()) && !(
                Constructible<T, optional<U> &       > ||
                Constructible<T, optional<U> &&      > ||
                Constructible<T, optional<U> const & > ||
                Constructible<T, optional<U> const &&> ||
                ConvertibleTo<optional<U> &,        T> ||
                ConvertibleTo<optional<U> &&,       T> ||
                ConvertibleTo<optional<U> const &,  T> ||
                ConvertibleTo<optional<U> const &&, T>
            )
    );

    CPP_def
    (
        template(typename U, typename T)
        concept OptionalShouldConvertAssign,
            OptionalShouldConvert<U, T> &&
            !(Assignable<T &, optional<U> &> ||
            Assignable<T &, optional<U> &&> ||
            Assignable<T &, optional<U> const &> ||
            Assignable<T &, optional<U> const &&>)
    );

    template<typename T>
    struct optional
      : detail::optional_adl::move_assign_layer<T>
    {
    private:
        using base_t = detail::optional_adl::move_assign_layer<T>;
    public:
        CPP_assert(Destructible<T>);
        static_assert(std::is_object<T>::value || std::is_lvalue_reference<T>::value, "");
        static_assert((bool) !Same<nullopt_t, uncvref_t<T>>, "");
        static_assert((bool) !Same<in_place_t, uncvref_t<T>>, "");
        using value_type = meta::_t<std::remove_cv<T>>;

        constexpr optional() noexcept
        {}
        constexpr optional(nullopt_t) noexcept
          : optional{}
        {}
        optional(optional const &) = default;
        optional(optional &&) = default;

        using base_t::base_t;

        CPP_template(typename E, typename... Args)(
            requires Constructible<T, std::initializer_list<E> &, Args...>)
        constexpr explicit optional(in_place_t, std::initializer_list<E> il, Args &&... args)
            noexcept(std::is_nothrow_constructible<T, std::initializer_list<E> &, Args...>::value)
          : base_t(in_place, il, static_cast<Args &&>(args)...)
        {}

        template<typename U = T>
        constexpr CPP_ctor(optional)(U &&v)(
            requires not defer::Same<detail::decay_t<U>, in_place_t> &&
                !defer::Same<detail::decay_t<U>, optional> &&
                defer::Constructible<T, U> &&
                defer::ConvertibleTo<U, T>)
          : base_t(in_place, static_cast<U &&>(v))
        {}
        template<typename U = T>
        explicit constexpr CPP_ctor(optional)(U &&v)(
            requires not defer::Same<detail::decay_t<U>, in_place_t> &&
                !defer::Same<detail::decay_t<U>, optional> &&
                defer::Constructible<T, U> &&
                !defer::ConvertibleTo<U, T>)
          : base_t(in_place, static_cast<U &&>(v))
        {}

        template<typename U>
        CPP_ctor(optional)(optional<U> const &that)(
            requires OptionalShouldConvert<U, T> &&
                Constructible<T, U const &> &&
                ConvertibleTo<U const &, T>)
        {
            if (that.has_value())
                base_t::construct_from(*that);
        }
        template<typename U>
        explicit CPP_ctor(optional)(optional<U> const &that)(
            requires OptionalShouldConvert<U, T> &&
                Constructible<T, U const &> &&
                !ConvertibleTo<U const &, T>)
        {
            if (that.has_value())
                base_t::construct_from(*that);
        }

        template<typename U>
        CPP_ctor(optional)(optional<U> &&that)(
            requires OptionalShouldConvert<U, T> &&
                Constructible<T, U> &&
                ConvertibleTo<U, T>)
        {
            if (that.has_value())
                base_t::construct_from(detail::move(*that));
        }
        template<typename U>
        explicit CPP_ctor(optional)(optional<U> &&that)(
            requires OptionalShouldConvert<U, T> &&
                Constructible<T, U> &&
                !ConvertibleTo<U, T>)
        {
            if (that.has_value())
                base_t::construct_from(detail::move(*that));
        }

        constexpr /*c++14*/
        optional &operator=(nullopt_t) noexcept
        {
            reset();
            return *this;
        }

        optional &operator=(optional const &) = default;
        optional &operator=(optional &&) = default;

        template<typename U = T>
        constexpr /*c++14*/
        auto operator=(U &&u)
            noexcept(std::is_nothrow_constructible<T, U>::value &&
                std::is_nothrow_assignable<T &, U>::value) ->
            CPP_ret(optional &)(
                requires not defer::Same<optional, detail::decay_t<U>> &&
                    !(defer::Satisfies<T, std::is_scalar> && defer::Same<T, detail::decay_t<U>>) &&
                    defer::Constructible<T, U> &&
                    defer::Assignable<T &, U>)
        {
            if (has_value())
                **this = static_cast<U &&>(u);
            else
                base_t::construct_from(static_cast<U &&>(u));
            return *this;
        }

        template<typename U>
        constexpr /*c++14*/
        auto operator=(optional<U> const &that) ->
            CPP_ret(optional &)(
                requires OptionalShouldConvertAssign<U, T> &&
                    Constructible<T, const U &> &&
                    Assignable<T &, const U &>)
        {
            base_t::assign_from(that);
            return *this;
        }

        template<typename U>
        constexpr /*c++14*/
        auto operator=(optional<U> &&that) ->
            CPP_ret(optional &)(
                requires OptionalShouldConvertAssign<U, T> &&
                    Constructible<T, U> &&
                    Assignable<T &, U>)
        {
            base_t::assign_from(std::move(that));
            return *this;
        }

        template<typename... Args>
        auto emplace(Args &&... args)
            noexcept(std::is_nothrow_constructible<T, Args...>::value) ->
            CPP_ret(T &)(
                requires Constructible<T, Args...>)
        {
            reset();
            return base_t::construct_from(static_cast<Args &&>(args)...);
        }
        template<typename E, typename... Args>
        auto emplace(std::initializer_list<E> il, Args &&... args)
            noexcept(std::is_nothrow_constructible<
                T, std::initializer_list<E> &, Args...>::value) ->
            CPP_ret(T &)(
                requires Constructible<T, std::initializer_list<E> &, Args &&...>)
        {
            reset();
            return base_t::construct_from(il, static_cast<Args &&>(args)...);
        }

        using base_t::swap;
        using base_t::operator->;
        using base_t::operator*;

        constexpr explicit operator bool() const noexcept
        {
            return has_value();
        }
        using base_t::has_value;

        constexpr T const &value() const &
        {
            return (has_value() || detail::throw_bad_optional_access()),
                **this;
        }
        constexpr /*c++14*/ T &value() &
        {
            return (has_value() || detail::throw_bad_optional_access()),
                **this;
        }
        constexpr T const &&value() const &&
        {
            return (has_value() || detail::throw_bad_optional_access()),
                detail::move(**this);
        }
        constexpr /*c++14*/ T &&value() &&
        {
            return (has_value() || detail::throw_bad_optional_access()),
                detail::move(**this);
        }

        template<typename U>
        constexpr auto value_or(U &&u) const & ->
            CPP_ret(T)(
                requires CopyConstructible<T> && ConvertibleTo<U, T>)
        {
            return has_value() ? **this : static_cast<T>((U &&)u);
        }
        template<typename U>
        constexpr auto value_or(U &&u) && ->
            CPP_ret(T)(
                requires MoveConstructible<T> && ConvertibleTo<U, T>)
        {
            return has_value() ? detail::move(**this) : static_cast<T>((U &&)u);
        }

        using base_t::reset;
    };

    /// \cond
    namespace detail
    {
        namespace optional_adl
        {
            constexpr bool convert_bool(bool b) noexcept
            {
                return b;
            }

            // Relational operators [optional.relops]
            template<typename T, typename U>
            constexpr auto operator==(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x == *y))) ->
                decltype(convert_bool(*x == *y))
            {
                return x.has_value() == y.has_value() &&
                    (!x || convert_bool(*x == *y));
            }
            template<typename T, typename U>
            constexpr auto operator!=(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x != *y))) ->
                decltype(convert_bool(*x != *y))
            {
                return x.has_value() != y.has_value() ||
                    (x && convert_bool(*x != *y));
            }
            template<typename T, typename U>
            constexpr auto operator<(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x < *y))) ->
                decltype(convert_bool(*x < *y))
            {
                return y && (!x || convert_bool(*x < *y));
            }
            template<typename T, typename U>
            constexpr auto operator>(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x > *y))) ->
                decltype(convert_bool(*x > *y))
            {
                return x && (!y || convert_bool(*x > *y));
            }
            template<typename T, typename U>
            constexpr auto operator<=(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x <= *y))) ->
                decltype(convert_bool(*x <= *y))
            {
                return !x || (y && convert_bool(*x <= *y));
            }
            template<typename T, typename U>
            constexpr auto operator>=(optional<T> const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(*x >= *y))) ->
                decltype(convert_bool(*x >= *y))
            {
                return !y || (x && convert_bool(*x >= *y));
            }

            // Comparisons with nullopt [optional.nullops]
            template<typename T>
            constexpr bool operator==(optional<T> const &x, nullopt_t) noexcept
            {
                return !x;
            }
            template<typename T>
            constexpr bool operator==(nullopt_t, optional<T> const &x) noexcept
            {
                return !x;
            }
            template<typename T>
            constexpr bool operator!=(optional<T> const &x, nullopt_t) noexcept
            {
                return !!x;
            }
            template<typename T>
            constexpr bool operator!=(nullopt_t, optional<T> const &x) noexcept
            {
                return !!x;
            }
            template<typename T>
            constexpr bool operator<(optional<T> const &, nullopt_t) noexcept
            {
                return false;
            }
            template<typename T>
            constexpr bool operator<(nullopt_t, optional<T> const &x) noexcept
            {
                return !!x;
            }
            template<typename T>
            constexpr bool operator>(optional<T> const &x, nullopt_t) noexcept
            {
                return !!x;
            }
            template<typename T>
            constexpr bool operator>(nullopt_t, optional<T> const &) noexcept
            {
                return false;
            }
            template<typename T>
            constexpr bool operator<=(optional<T> const &x, nullopt_t) noexcept
            {
                return !x;
            }
            template<typename T>
            constexpr bool operator<=(nullopt_t, optional<T> const &) noexcept
            {
                return true;
            }
            template<typename T>
            constexpr bool operator>=(optional<T> const &, nullopt_t) noexcept
            {
                return true;
            }
            template<typename T>
            constexpr bool operator>=(nullopt_t, optional<T> const &x) noexcept
            {
                return !x;
            }

            // Comparisons with T [optional.comp_with_t]
            template<typename T, typename U>
            constexpr auto operator==(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x == y))) ->
                decltype(convert_bool(*x == y))
            {
                return x && convert_bool(*x == y);
            }
            template<typename T, typename U>
            constexpr auto operator==(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x == *y))) ->
                decltype(convert_bool(x == *y))
            {
                return y && convert_bool(x == *y);
            }
            template<typename T, typename U>
            constexpr auto operator!=(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x != y))) ->
                decltype(convert_bool(*x != y))
            {
                return !x || convert_bool(*x != y);
            }
            template<typename T, typename U>
            constexpr auto operator!=(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x != *y))) ->
                decltype(convert_bool(x != *y))
            {
                return !y || convert_bool(x != *y);
            }
            template<typename T, typename U>
            constexpr auto operator<(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x < y))) ->
                decltype(convert_bool(*x < y))
            {
                return !x || convert_bool(*x < y);
            }
            template<typename T, typename U>
            constexpr auto operator<(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x < *y))) ->
                decltype(convert_bool(x < *y))
            {
                return y && convert_bool(x < *y);
            }
            template<typename T, typename U>
            constexpr auto operator>(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x > y))) ->
                decltype(convert_bool(*x > y))
            {
                return x && convert_bool(*x > y);
            }
            template<typename T, typename U>
            constexpr auto operator>(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x > *y))) ->
                decltype(convert_bool(x > *y))
            {
                return !y || convert_bool(x > *y);
            }
            template<typename T, typename U>
            constexpr auto operator<=(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x <= y))) ->
                decltype(convert_bool(*x <= y))
            {
                return !x || convert_bool(*x <= y);
            }
            template<typename T, typename U>
            constexpr auto operator<=(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x <= *y))) ->
                decltype(convert_bool(x <= *y))
            {
                return y && convert_bool(x <= *y);
            }
            template<typename T, typename U>
            constexpr auto operator>=(optional<T> const &x, U const &y)
                noexcept(noexcept(convert_bool(*x >= y))) ->
                decltype(convert_bool(*x >= y))
            {
                return x && convert_bool(*x >= y);
            }
            template<typename T, typename U>
            constexpr auto operator>=(T const &x, optional<U> const &y)
                noexcept(noexcept(convert_bool(x >= *y))) ->
                decltype(convert_bool(x >= *y))
            {
                return !y || convert_bool(x >= *y);
            }

            template<typename T>
            auto CPP_auto_fun(swap)(optional<T> &x, optional<T> &y)
            (
                return x.swap(y)
            )
        } // namespace optional_adl
    } // namespace detail
    /// \endcond

    template<typename T>
    constexpr auto CPP_auto_fun(make_optional)(T &&t)
    (
        return optional<detail::decay_t<T>>{static_cast<T &&>(t)}
    )
    template<typename T, typename... Args>
    constexpr auto CPP_auto_fun(make_optional)(Args &&... args)
    (
        return optional<T>{in_place, static_cast<Args &&>(args)...}
    )
    template<typename T, typename U, typename... Args>
    constexpr auto CPP_auto_fun(make_optional)(std::initializer_list<U> il, Args &&... args)
    (
        return optional<T>{in_place, il, static_cast<Args &&>(args)...}
    )

    /// \cond
    namespace detail
    {
        template<typename T, typename Tag = void, bool Enable = true>
        struct non_propagating_cache
          : optional<T>
        {
            non_propagating_cache() = default;
            constexpr non_propagating_cache(nullopt_t) noexcept
            {}
            constexpr
            non_propagating_cache(non_propagating_cache const &) noexcept
              : optional<T>{}
            {}
            constexpr /*c++14*/
            non_propagating_cache(non_propagating_cache &&that) noexcept
              : optional<T>{}
            {
                that.optional<T>::reset();
            }
            constexpr /*c++14*/
            non_propagating_cache &operator=(non_propagating_cache const &) noexcept
            {
                optional<T>::reset();
                return *this;
            }
            constexpr /*c++14*/
            non_propagating_cache &operator=(non_propagating_cache &&that) noexcept
            {
                that.optional<T>::reset();
                optional<T>::reset();
                return *this;
            }
            using optional<T>::operator=;
        };

        template<typename T, typename Tag>
        struct non_propagating_cache<T, Tag, false>
        {};
    } // namespace detail
    /// \endcond
} // namespace ranges

#endif
