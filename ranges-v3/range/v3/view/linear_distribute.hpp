/// \file
// Range v3 library
//
//  Copyright Casey Carter 2017
//  Copyright Gonzalo Brito Gadeschi 2017
//
//  Use, modification and distribution is subject to the
//  Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// Project home: https://github.com/ericniebler/range-v3
//

#ifndef RANGES_V3_VIEW_LINEAR_DISTRIBUTE_HPP
#define RANGES_V3_VIEW_LINEAR_DISTRIBUTE_HPP

#include <type_traits>
#include <range/v3/range_fwd.hpp>
#include <range/v3/view/facade.hpp>
#include <range/v3/utility/static_const.hpp>
#include <range/v3/iterator/default_sentinel.hpp>

namespace ranges
{
    namespace view
    {
        /// \addtogroup group-views
        /// @{

        template<typename T>
        struct linear_distribute_view
          : view_facade<linear_distribute_view<T>, finite>
        {
            CPP_assert(std::is_arithmetic<T>());
        private:
            friend range_access;

            T from_, to_;
            std::ptrdiff_t n_;

            constexpr T read() const noexcept
            {
                return from_;
            }
            constexpr bool equal(default_sentinel_t) const noexcept
            {
                return n_ == 0;
            }
            constexpr /*c++14*/
            bool equal(linear_distribute_view const& other) const noexcept
            {
                RANGES_DIAGNOSTIC_PUSH
                RANGES_DIAGNOSTIC_IGNORE_FLOAT_EQUAL
                RANGES_EXPECT(from_ == other.from_ && to_ == other.to_);
                RANGES_DIAGNOSTIC_POP
                return n_ == other.n_;
            }
            constexpr /*c++14*/ void next() noexcept
            {
                RANGES_EXPECT(n_ > 0);
                --n_;
                if(n_ == 0)
                {
                    from_ = to_;
                }
                else
                {
                    from_ += (to_ - from_) / T(n_);
                }
            }
        public:
            constexpr /*c++14*/ linear_distribute_view() = default;
            constexpr /*c++14*/ linear_distribute_view(T from, T to__, std::ptrdiff_t n) noexcept
              : from_(from), to_(to__), n_(n)
            {
                RANGES_EXPECT(n_ > 0);
                RANGES_EXPECT(to_ >= from_);
            }
            constexpr std::size_t size() const noexcept
            {
                return static_cast<std::size_t>(n_);
            }
        };

        /// Distributes `n` values linearly in the closed interval [`from`, `to`].
        ///
        /// \pre `from <= to && n > 0`
        ///
        /// If `from == to`, returns n-times `to`.
        /// If `n == 1` returns `to`.
        struct linear_distribute_fn
        {
            template<typename T>
            constexpr auto CPP_fun(operator())(T from, T to, std::ptrdiff_t n) (const
                requires std::is_arithmetic<T>::value)
            {
                return linear_distribute_view<T>{from, to, n};
            }
        };

        /// \relates linear_distribute_fn
        /// \ingroup group-views
        RANGES_INLINE_VARIABLE(linear_distribute_fn, linear_distribute)
    }
}

#endif

