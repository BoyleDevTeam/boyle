/**
 * @file utils.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "spdlog/spdlog.h"

#include "math/type_traits.hpp"
#include "math/vec2.hpp"
#include "math/vec3.hpp"

namespace tiny_pnc {
namespace math {

constexpr double kEpsilon{1e-8};

template <typename T>
[[using gnu: const, always_inline, leaf]]
constexpr inline T pow(T x, std::size_t n) noexcept {
    static_assert(std::is_arithmetic_v<T>, "The loaded type must support arithmetic operations.");
    T result{1.0};
    for (std::size_t i{0}; i < n; ++i) {
        result *= x;
    }
    return result;
}

template <typename T>
[[using gnu: const, always_inline, leaf]]
constexpr inline bool inRange(T value, T start, T end) noexcept {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    if (start > end) {
        start -= end;
        end += start;
        start = end - start;
    }
    return value > start && value < end;
}

template <typename T>
[[using gnu: const, flatten, leaf]] [[nodiscard]]
inline std::vector<T> linspace(T start, T end, std::size_t num, bool endpoint = true) noexcept {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    if (num == 0) {
        return std::vector<T>(num);
    } else if (num == 1) {
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

template <typename T, typename U = double>
[[using gnu: const, always_inline, leaf, hot]] [[nodiscard]]
constexpr inline T lerp(T start, T end, U ratio) noexcept {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    return (1 - ratio) * start + ratio * end;
}

template <
    typename ForwardIt, typename T = typename std::iterator_traits<ForwardIt>::value_type,
    typename U = double, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
[[using gnu: pure, flatten, leaf]]
inline bool hasDuplicates(ForwardIt first, ForwardIt last, U tol = kEpsilon) {
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    if (last - first < 2) {
        std::string error_msg = std::format(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    if (std::is_sorted(first, last)) {
        for (ForwardIt it{first}; it != last - 1; ++it) {
            if (std::abs(*(it + 1) - *it) < tol) {
                return true;
            }
        }
    } else {
        std::vector<T> arr{first, last};
        std::sort(arr.begin(), arr.end());
        for (typename std::vector<T>::const_iterator it{arr.cbegin()}; it != arr.cend() - 1; ++it) {
            if (std::abs(*(it + 1) - *it) < tol) {
                return true;
            }
        }
    }
    return false;
}

template <
    typename ForwardIt, typename T = typename std::iterator_traits<ForwardIt>::value_type,
    typename U = typename T::value_type,
    std::enable_if_t<tiny_pnc::math::isVecArithmeticV<T>, bool> = true>
[[using gnu: pure, flatten, leaf]]
inline bool hasDuplicates(ForwardIt first, ForwardIt last, U tol = kEpsilon) {
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    if (last - first < 2) {
        std::string error_msg = std::format(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    for (ForwardIt it{first}; it != last - 1; ++it) {
        if (std::hypot((*(it + 1) - *it)) < tol) {
            return true;
        }
    }
    return false;
}

template <
    typename ForwardIt, typename T = typename std::iterator_traits<ForwardIt>::value_type,
    typename U = typename T::value_type,
    std::enable_if_t<std::is_same_v<T, Vec2<Vec2Mode::XY, U>>, bool> = true>
[[using gnu: pure, flatten, leaf]]
inline ForwardIt closetUpperElement(ForwardIt first, ForwardIt last, T element) {
    if (last - first < 2) {
        std::string error_msg = std::format(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    const std::size_t size = last - first;
    std::size_t pos;
    U min_distance = std::numeric_limits<U>::max();
    for (ForwardIt it{first}; it != last; ++it) {
        const U distance = element.distanceTo(*it);
        if (distance < min_distance) {
            pos = it - first;
            min_distance = distance;
        }
    }
    T diff, r;
    U inner_prod;
    if (pos == 0) {
        diff = *(first + 1) - *first;
        r = element - *first;
        inner_prod = diff.dot(r);
        if (inner_prod < 0.0) {
            return first;
        } else {
            return first + 1;
        }
    } else if (pos == size - 1) {
        diff = *(last - 2) - *(last - 1);
        r = element - *(last - 1);
        inner_prod = diff.dot(r);
        if (inner_prod < 0.0) {
            return last;
        } else {
            return last - 1;
        }
    } else {
        ForwardIt it = first + pos;
        diff = *(it + 1) - *it;
        r = element - *it;
        inner_prod = diff.dot(r);
        if (inner_prod < 0.0) {
            return it;
        } else {
            return it + 1;
        }
    }
}

template <Vec2Mode VM, typename T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline std::vector<Vec2<VM, T>> stack(const std::vector<T>& xs, const std::vector<T>& ys) {
    if (xs.size() != ys.size()) {
        std::string error_msg = std::format(
            "Invalid arguments detected! xs, ys must share the same size: xs.size() = {0:d} while "
            "ys.size() = {1:d}",
            xs.size(), ys.size()
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    const std::size_t size{xs.size()};
    std::vector<Vec2<VM, T>> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i]);
    }
    return result;
}

template <Vec3Mode VM, typename T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline std::vector<Vec3<VM, T>> stack(
    const std::vector<T>& xs, const std::vector<T>& ys, const std::vector<T>& zs
) {
    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        std::string error_msg = std::format(
            "Invalid arguments detected! xs, ys must share the same size: xs.size() = {0:d} while "
            "ys.size() = {1:d}",
            xs.size(), ys.size()
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    const std::size_t size{xs.size()};
    std::vector<Vec3<VM, T>> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i], zs[i]);
    }
    return result;
}

} // namespace math
} // namespace tiny_pnc
