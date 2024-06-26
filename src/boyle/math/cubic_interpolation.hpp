/**
 * @file cubic_interpolation.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <concepts>

#include "boyle/math/concepts.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto cuberp(T start, T end, T ddstart, T ddend, U ratio, U scale = 1.0) noexcept
    -> T {
    constexpr std::array<double, 2> kFactors{-(1.0 / 6.0), 1.0 / 6.0};
    const U ratio3{ratio * ratio * ratio};
    const U resi_ratio{1.0 - ratio};
    const U resi_ratio3{resi_ratio * resi_ratio * resi_ratio};
    return lerp(start, end, ratio) +
           ((resi_ratio * kFactors[0] + resi_ratio3 * kFactors[1]) * ddstart +
            (ratio * kFactors[0] + ratio3 * kFactors[1]) * ddend) *
               scale * scale;
}

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, leaf, hot]] [[nodiscard]]
inline constexpr auto cuberpd(T start, T end, T ddstart, T ddend, U ratio, U scale = 1.0) noexcept
    -> T {
    constexpr std::array<U, 2> kFactors{-(1.0 / 6.0), 0.5};
    const U ratio2{ratio * ratio};
    const U resi_ratio{1.0 - ratio};
    const U resi_ratio2{resi_ratio * resi_ratio};
    return (end - start) / scale + (-(kFactors[0] + resi_ratio2 * kFactors[1]) * ddstart +
                                    (kFactors[0] + ratio2 * kFactors[1]) * ddend) *
                                       scale;
}

template <std::floating_point T = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto cuberpCoeffs(T ratio, T scale = 1.0) noexcept -> std::array<T, 4> {
    constexpr std::array<T, 2> kFactors{-(1.0 / 6.0), 1.0 / 6.0};
    const T ratio3{ratio * ratio * ratio};
    const T resi_ratio{1.0 - ratio};
    const T resi_ratio3{resi_ratio * resi_ratio * resi_ratio};
    const T scale2{scale * scale};
    return std::array<T, 4>{
        resi_ratio, ratio, (resi_ratio * kFactors[0] + resi_ratio3 * kFactors[1]) * scale2,
        (ratio * kFactors[0] + ratio3 * kFactors[1]) * scale2
    };
}

template <std::floating_point T = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto cuberpdCoeffs(T ratio, T scale = 1.0) noexcept -> std::array<T, 4> {
    constexpr std::array<T, 2> kFactors{-(1.0 / 6.0), 0.5};
    const T ratio2{ratio * ratio};
    const T resi_ratio{1.0 - ratio};
    const T resi_ratio2{resi_ratio * resi_ratio};
    const T reci_scale{1.0 / scale};
    return std::array<T, 4>{
        -reci_scale, reci_scale, -(kFactors[0] + resi_ratio2 * kFactors[1]) * scale,
        (kFactors[0] + ratio2 * kFactors[1]) * scale
    };
}

} // namespace boyle::math
