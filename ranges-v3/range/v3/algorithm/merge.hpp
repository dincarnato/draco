/// \file
// Range v3 library
//
//  Copyright Eric Niebler 2014-present
//
//  Use, modification and distribution is subject to the
//  Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// Project home: https://github.com/ericniebler/range-v3
//
// Copyright (c) 2009 Alexander Stepanov and Paul McJones
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without
// fee, provided that the above copyright notice appear in all copies
// and that both that copyright notice and this permission notice
// appear in supporting documentation. The authors make no
// representations about the suitability of this software for any
// purpose. It is provided "as is" without express or implied
// warranty.
//
// Algorithms from
// Elements of Programming
// by Alexander Stepanov and Paul McJones
// Addison-Wesley Professional, 2009
#ifndef RANGES_V3_ALGORITHM_MERGE_HPP
#define RANGES_V3_ALGORITHM_MERGE_HPP

#include <tuple>
#include <range/v3/range_fwd.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/dangling.hpp>
#include <range/v3/range/traits.hpp>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/result_types.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/functional/invoke.hpp>
#include <range/v3/iterator/traits.hpp>
#include <range/v3/utility/static_const.hpp>

namespace ranges
{
    /// \addtogroup group-algorithms
    /// @{
    template<typename I0, typename I1, typename O>
    using merge_result = detail::in1_in2_out_result<I0, I1, O>;

    struct merge_fn
    {
        template<typename I0, typename S0, typename I1, typename S1, typename O,
            typename C = less, typename P0 = identity, typename P1 = identity>
        auto operator()(I0 begin0, S0 end0, I1 begin1, S1 end1, O out, C pred = C{},
                P0 proj0 = P0{}, P1 proj1 = P1{}) const ->
            CPP_ret(merge_result<I0, I1, O>)(
                requires Sentinel<S0, I0> && Sentinel<S1, I1> &&
                    Mergeable<I0, I1, O, C, P0, P1>)
        {
            for(; begin0 != end0 && begin1 != end1; ++out)
            {
                if(invoke(pred, invoke(proj1, *begin1), invoke(proj0, *begin0)))
                {
                    *out = *begin1;
                    ++begin1;
                }
                else
                {
                    *out = *begin0;
                    ++begin0;
                }
            }
            auto t0 = copy(begin0, end0, out);
            auto t1 = copy(begin1, end1, t0.out);
            return {t0.in, t1.in, t1.out};
        }

        template<typename Rng0, typename Rng1, typename O, typename C = less,
            typename P0 = identity, typename P1 = identity>
        auto operator()(Rng0 &&rng0, Rng1 &&rng1, O out, C pred = C{}, P0 proj0 = P0{},
                P1 proj1 = P1{}) const ->
            CPP_ret(merge_result<safe_iterator_t<Rng0>, safe_iterator_t<Rng1>, O>)(
                requires Range<Rng0> && Range<Rng1> &&
                    Mergeable<iterator_t<Rng0>, iterator_t<Rng1>, O, C, P0, P1>)
        {
            return (*this)(begin(rng0), end(rng0), begin(rng1), end(rng1), std::move(out),
                std::move(pred), std::move(proj0), std::move(proj1));
        }
    };

    /// \sa `merge_fn`
    /// \ingroup group-algorithms
    RANGES_INLINE_VARIABLE(merge_fn, merge)
    /// @}
} // namespace ranges

#endif // include guard
