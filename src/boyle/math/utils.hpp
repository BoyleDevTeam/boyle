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
#include <ranges>
#include <utility>
#include <vector>

#include "boyle/math/concepts.hpp"

namespace boyle::math {

inline constexpr double kEpsilon{1e-8};

struct periodic_tag final {};

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
    Arithmetic auto value, Arithmetic auto start, Arithmetic auto end, double tol = 1E-8
) noexcept -> bool {
    const double c_value = value;
    return (c_value - start) * (c_value - end) < -tol;
}

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, leaf, hot]] [[nodiscard]]
inline constexpr auto lerp(T start, T end, U ratio) noexcept -> T {
    return (1 - ratio) * start + ratio * end;
}

template <GeneralArithmetic T>
[[using gnu: const, always_inline]] [[nodiscard]]
inline auto linspace(T start, T end, std::size_t num, bool endpoint = true) noexcept
    -> std::vector<T> {
    if (num == 0) [[unlikely]] {
        return std::vector<T>{};
    }
    if (num == 1) [[unlikely]] {
        return std::vector<T>{start};
    }
    std::vector<T> result(num);
    T step = endpoint ? (end - start) / static_cast<double>(num - 1)
                      : (end - start) / static_cast<double>(num);
    std::ranges::generate(result, [t = start - step, h = step]() mutable constexpr noexcept -> T {
        return t += h;
    });
    return result;
}

[[using gnu: pure, always_inline]]
inline auto hasDuplicates(std::ranges::forward_range auto&& range, double tol = 1E-8) noexcept
    -> bool
    requires Arithmetic<std::ranges::range_value_t<decltype(range)>>
{
    if (range.size() < 2) [[unlikely]] {
        return false;
    }
    std::vector<std::ranges::range_value_t<decltype(range)>> sorted{range.begin(), range.end()};
    std::sort(sorted.begin(), sorted.end());
    for (auto pre{sorted.cbegin()}, cur{std::next(sorted.cbegin())}; cur != sorted.cend();
         ++pre, ++cur) {
        if (std::abs(*cur - *pre) < tol) {
            return true;
        }
    }
    return false;
}

[[using gnu: pure, always_inline]]
inline auto hasDuplicates(std::ranges::forward_range auto&& range, double tol = 1E-8) noexcept
    -> bool
    requires VecArithmetic<std::ranges::range_value_t<decltype(range)>>
{
    if (range.size() < 2) [[unlikely]] {
        return false;
    }
    for (auto pre{range.begin()}, cur{std::next(range.begin())}; cur != range.end(); ++pre, ++cur) {
        if (std::abs(*cur - *pre) < tol) {
            return true;
        }
    }
    return false;
}

[[using gnu: pure, flatten, leaf, hot]]
inline auto nearestUpperElement(
    std::ranges::forward_range auto&& range, std::ranges::range_value_t<decltype(range)> element,
    double tol = 1E-8
) noexcept -> std::ranges::iterator_t<decltype(range)>
    requires Arithmetic<std::ranges::range_value_t<decltype(range)>>
{
    if (range.size() < 2) [[unlikely]] {
        return element < range.front() ? range.begin() : range.end();
    }
    if (std::abs(element - range.front()) < tol) {
        return std::next(range.begin());
    }
    if (std::abs(element - *(std::prev(range.end()))) < tol) {
        return std::prev(range.end());
    }
    return std::ranges::upper_bound(
        std::forward<decltype(range)>(range), element,
        [tol](
            const std::ranges::range_value_t<decltype(range)>& lhs,
            const std::ranges::range_value_t<decltype(range)>& rhs
        ) noexcept -> bool { return lhs - rhs < -tol; }
    );
}

[[using gnu: pure, flatten, leaf, hot]]
inline auto nearestUpperElement(
    std::ranges::forward_range auto&& range, std::ranges::range_value_t<decltype(range)> element,
    double tol = 1E-8
) noexcept -> std::ranges::iterator_t<decltype(range)>
    requires VecArithmetic<std::ranges::range_value_t<decltype(range)>>
{
    if (range.size() < 2) [[unlikely]] {
        return range.begin();
    }
    std::size_t pos{0};
    double min_euclidean = std::numeric_limits<double>::max();
    for (auto it{range.begin()}; it != range.end(); ++it) {
        const double euclidean = element.euclideanTo(*it);
        if (euclidean < min_euclidean) {
            pos = it - range.begin();
            min_euclidean = euclidean;
        }
    }
    std::ranges::range_value_t<decltype(range)> diff, r;
    double inner_prod{0.0};
    if (pos == 0) {
        diff = *(std::next(range.begin())) - range.front();
        r = element - range.front();
        inner_prod = diff.dot(r);
        return inner_prod < -tol ? range.begin() : std::next(range.begin());
    }
    if (pos == range.size() - 1) {
        diff = *(std::prev(range.end(), 2)) - *(std::prev(range.end()));
        r = element - *(std::prev(range.end()));
        inner_prod = diff.dot(r);
        return inner_prod < -tol ? range.end() : std::prev(range.end());
    }
    std::ranges::iterator_t<decltype(range)> it = range.begin() + pos;
    diff = *(std::next(it)) - *it;
    r = element - *it;
    inner_prod = diff.dot(r);
    return inner_prod < -tol ? it : std::next(it);
}

} // namespace boyle::math
