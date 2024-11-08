/// \file
// Range v3 library
//
//  Copyright Eric Niebler 2013-present
//
//  Use, modification and distribution is subject to the
//  Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// Project home: https://github.com/ericniebler/range-v3
//

#ifndef RANGES_V3_VIEW_SINGLE_HPP
#define RANGES_V3_VIEW_SINGLE_HPP

#include <utility>
#include <type_traits>
#include <range/v3/range_fwd.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/traits.hpp>
#include <range/v3/view/facade.hpp>
#include <range/v3/iterator/concepts.hpp>
#include <range/v3/iterator/traits.hpp>
#include <range/v3/utility/optional.hpp>
#include <range/v3/utility/semiregular.hpp>
#include <range/v3/utility/static_const.hpp>

namespace ranges
{
    /// \addtogroup group-views
    /// @{
    template<typename T>
    struct single_view
      : view_interface<single_view<T>, (cardinality)1> {
    private:
        CPP_assert(CopyConstructible<T>);
        static_assert(std::is_object<T>::value,
            "The template parameter of single_view must be an object type");
        semiregular_t<T> value_;
        template<typename... Args>
        constexpr single_view(in_place_t, std::true_type, Args &&...args)
          : value_{static_cast<Args &&>(args)...}
        {}
        template<typename... Args>
        constexpr single_view(in_place_t, std::false_type, Args &&...args)
          : value_{in_place, static_cast<Args &&>(args)...}
        {}
    public:
        single_view() = default;
        constexpr explicit single_view(T const &t)
          : value_(t)
        {}
        constexpr explicit single_view(T &&t)
          : value_(std::move(t))
        {}
        CPP_template(class... Args)(
            requires Constructible<T, Args...>)
        constexpr single_view(in_place_t, Args&&... args)
          : single_view{
                in_place,
                meta::bool_<(bool) Semiregular<T>>{},
                static_cast<Args &&>(args)...}
        {}
        constexpr /*c++14*/ T *begin() noexcept
        {
            return data();
        }
        constexpr T const *begin() const noexcept
        {
            return data();
        }
        constexpr /*c++14*/ T *end() noexcept
        {
            return data() + 1;
        }
        constexpr T const *end() const noexcept
        {
            return data() + 1;
        }
        constexpr static std::size_t size() noexcept
        {
            return 1u;
        }
        constexpr /*c++14*/ T *data() noexcept
        {
            return std::addressof(static_cast<T &>(value_));
        }
        constexpr T const *data() const noexcept
        {
            return std::addressof(static_cast<T const &>(value_));
        }
    };

#if RANGES_CXX_DEDUCTION_GUIDES >= RANGES_CXX_DEDUCTION_GUIDES_17
    template<class T>
    explicit single_view(T &&) -> single_view<detail::decay_t<T>>;
#endif

    namespace view
    {
        struct single_fn
        {
            template<typename Val>
            auto operator()(Val value) const ->
                CPP_ret(single_view<Val>)(
                    requires CopyConstructible<Val>)
            {
                return single_view<Val>{std::move(value)};
            }
        };

        /// \relates single_fn
        /// \ingroup group-views
        RANGES_INLINE_VARIABLE(single_fn, single)
    }
    /// @}
}

#include <range/v3/detail/satisfy_boost_range.hpp>
RANGES_SATISFY_BOOST_RANGE(::ranges::single_view)

#endif
