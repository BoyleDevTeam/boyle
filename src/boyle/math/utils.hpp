/**
 * @file utils.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <iterator>
#include <vector>

#include "boyle/common/utils/logging.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/duplet.hpp"
#include "boyle/math/triplet.hpp"

namespace boyle::math {

inline constexpr double kEpsilon{1e-8};

[[using gnu: const, always_inline, leaf]]
inline constexpr auto pow(Arithmetic auto x, std::size_t n) noexcept -> decltype(x) {
    decltype(x) result = 1.0;
    for (; n > 0; --n) {
        result *= x;
    }
    return result;
}

[[using gnu: const, always_inline, leaf]]
inline constexpr auto inRange(
    Arithmetic auto value, Arithmetic auto start, Arithmetic auto end, double tol = kEpsilon
) noexcept -> bool {
    const double c_value = value;
    return (c_value - start) * (c_value - end) < -tol;
}

template <GeneralArithmetic T>
[[using gnu: const, always_inline, leaf]] [[nodiscard]]
inline auto linspace(T start, T end, std::size_t num, bool endpoint = true) noexcept
    -> std::vector<T> {
    if (num == 0) {
        return std::vector<T>{};
    }
    if (num == 1) {
        return std::vector<T>(1, start);
    }
    std::vector<T> result;
    result.reserve(num);
    T step;
    if (endpoint) {
        step = (end - start) / static_cast<double>(num - 1);
    } else {
        step = (end - start) / static_cast<double>(num);
        end -= step;
    }
    for (T x{start}; num; x += step, --num) {
        result.push_back(x);
    }
    return result;
}

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, leaf, hot]] [[nodiscard]]
inline constexpr auto lerp(T start, T end, U ratio) noexcept -> T {
    return (1 - ratio) * start + ratio * end;
}

template <std::forward_iterator ForwardIt>
[[using gnu: const, always_inline, hot]]
inline auto hasDuplicates(
    ForwardIt first, std::sentinel_for<ForwardIt> auto last, double tol = kEpsilon
) noexcept -> bool
    requires Arithmetic<typename std::iterator_traits<ForwardIt>::value_type>
{
    if (last - first < 2) {
        BOYLE_LOG_WARN(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        return false;
    }
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    std::vector<T> sorted{first, last};
    std::sort(sorted.begin(), sorted.end());
    for (typename std::vector<T>::const_iterator it{sorted.cbegin()}; it != sorted.cend() - 1;
         ++it) {
        if (std::abs(*(it + 1) - *it) < tol) {
            return true;
        }
    }
    return false;
}

template <std::forward_iterator ForwardIt>
[[using gnu: const, always_inline, hot]]
inline auto hasDuplicates(
    ForwardIt first, std::sentinel_for<ForwardIt> auto last, double tol = kEpsilon
) noexcept -> bool
    requires VecArithmetic<typename std::iterator_traits<ForwardIt>::value_type>
{
    if (last - first < 2) {
        BOYLE_LOG_WARN(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        return false;
    }
    for (ForwardIt it{first}; it != last - 1; ++it) {
        if (std::hypot((*(it + 1) - *it)) < tol) {
            return true;
        }
    }
    return false;
}

template <std::forward_iterator ForwardIt>
[[using gnu: pure, flatten, leaf, hot]]
inline auto nearestUpperElement(
    ForwardIt first, std::sentinel_for<ForwardIt> auto last,
    typename std::iterator_traits<ForwardIt>::value_type element, double tol = kEpsilon
) noexcept -> decltype(first)
    requires Arithmetic<typename std::iterator_traits<ForwardIt>::value_type>
{
    if (last - first < 2) {
        BOYLE_LOG_WARN(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        return element < *first ? first : last;
    }
    if (std::abs(element - *first) < tol) {
        return first + 1;
    }
    if (std::abs(element - *(last - 1)) < tol) {
        return last - 1;
    }
    return std::ranges::upper_bound(
        first, last, element,
        [tol](
            typename std::iterator_traits<ForwardIt>::value_type lhs,
            typename std::iterator_traits<ForwardIt>::value_type rhs
        ) -> bool { return lhs - rhs < -tol; }
    );
}

template <std::forward_iterator ForwardIt>
[[using gnu: pure, flatten, leaf, hot]]
inline auto nearestUpperElement(
    ForwardIt first, std::sentinel_for<ForwardIt> auto last,
    typename std::iterator_traits<ForwardIt>::value_type element, double tol = kEpsilon
) noexcept -> ForwardIt
    requires VecArithmetic<typename std::iterator_traits<ForwardIt>::value_type>
{
    if (last - first < 2) {
        BOYLE_LOG_ERROR(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        return ForwardIt{};
    }
    const std::size_t size = last - first;
    std::size_t pos{0};
    double min_euclidean = std::numeric_limits<double>::max();
    for (ForwardIt it{first}; it != last; ++it) {
        const double euclidean = element.euclideanTo(*it);
        if (euclidean < min_euclidean) {
            pos = it - first;
            min_euclidean = euclidean;
        }
    }
    typename std::iterator_traits<ForwardIt>::value_type diff, r;
    double inner_prod{0.0};
    if (pos == 0) {
        diff = *(first + 1) - *first;
        r = element - *first;
        inner_prod = diff.dot(r);
        return inner_prod < -tol ? first : first + 1;
    }
    if (pos == size - 1) {
        diff = *(last - 2) - *(last - 1);
        r = element - *(last - 1);
        inner_prod = diff.dot(r);
        return inner_prod < -tol ? last : last - 1;
    }
    ForwardIt it = first + pos;
    diff = *(it + 1) - *it;
    r = element - *it;
    inner_prod = diff.dot(r);
    return inner_prod < -tol ? it : it + 1;
}

template <Dupletic T>
[[using gnu: pure, always_inline, leaf]] [[nodiscard]]
inline auto squeeze(
    const std::vector<DupletElementT<0, T>>& xs, const std::vector<DupletElementT<1, T>>& ys
) noexcept -> std::vector<T> {
    if (xs.size() != ys.size()) {
        BOYLE_LOG_ERROR(
            "Invalid arguments detected! xs, ys must share the same size: xs.size() = {0:d} while "
            "ys.size() = {1:d}.",
            xs.size(), ys.size()
        );
        return std::vector<T>{};
    }
    const std::size_t size{xs.size()};
    std::vector<T> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i]);
    }
    return result;
}

template <Tripletic T>
[[using gnu: pure, always_inline, leaf]] [[nodiscard]]
inline auto squeeze(
    const std::vector<TripletElementT<0, T>>& xs, const std::vector<TripletElementT<1, T>>& ys,
    const std::vector<TripletElementT<2, T>>& zs
) noexcept -> std::vector<T> {
    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        BOYLE_LOG_ERROR(
            "Invalid arguments detected! xs, ys, zs must share the same size: xs.size() = {0:d} "
            "while ys.size() = {1:d} while zs.size() = {2:d}.",
            xs.size(), ys.size(), zs.size()
        );
        return std::vector<T>{};
    }
    const std::size_t size{xs.size()};
    std::vector<T> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i], zs[i]);
    }
    return result;
}

DECLARE_DUPLET(SlDuplet, s, l);

using SlDupletf = SlDuplet<float>;
using SlDupletd = SlDuplet<double>;

DECLARE_DUPLET(StDuplet, s, t);

using StDupletf = StDuplet<float>;
using StDupletd = StDuplet<double>;

DECLARE_DUPLET(SxDuplet, s, x);

using SxDupletf = SxDuplet<float>;
using SxDupletd = SxDuplet<double>;

DECLARE_DUPLET(SyDuplet, s, y);

using SyDupletf = SyDuplet<float>;
using SyDupletd = SyDuplet<double>;

DECLARE_TRIPLET(SlzTriplet, s, l, z);

using SlzTripletf = SlzTriplet<float>;
using SlzTripletd = SlzTriplet<double>;

} // namespace boyle::math
