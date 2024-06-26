/**
 * @file quintic_interpolation.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <concepts>

#include "boyle/math/concepts.hpp"
#include "boyle/math/cubic_interpolation.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto quinerp(
    T start, T end, T ddstart, T ddend, T d4start, T d4end, U ratio, U scale = 1.0
) noexcept -> T {
    constexpr std::array<U, 3> kFactors{7.0 / 360.0, -(1.0 / 36.0), 1.0 / 120.0};
    const U ratio3{ratio * ratio * ratio};
    const U ratio5{ratio3 * ratio * ratio};
    const U resi_ratio{1.0 - ratio};
    const U resi_ratio3{resi_ratio * resi_ratio * resi_ratio};
    const U resi_ratio5{resi_ratio3 * resi_ratio * resi_ratio};
    return cuberp(start, end, ddstart, ddend, ratio, scale) +
           ((resi_ratio * kFactors[0] + resi_ratio3 * kFactors[1] + resi_ratio5 * kFactors[2]) *
                d4start +
            (ratio * kFactors[0] + ratio3 * kFactors[1] + ratio5 * kFactors[2]) * d4end) *
               scale * scale * scale * scale;
}

template <GeneralArithmetic T, std::floating_point U = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto quinerpd(
    T start, T end, T ddstart, T ddend, T d4start, T d4end, U ratio, U scale = 1.0
) noexcept -> T {
    constexpr std::array<U, 3> kFactors{7.0 / 360.0, -(1.0 / 12.0), 1.0 / 24.0};
    const U ratio2{ratio * ratio};
    const U ratio4{ratio2 * ratio2};
    const U resi_ratio{1.0 - ratio};
    const U resi_ratio2{resi_ratio * resi_ratio};
    const U resi_ratio4{resi_ratio2 * resi_ratio2};
    return cuberpd(start, end, ddstart, ddend, ratio, scale) +
           (-(kFactors[0] + resi_ratio2 * kFactors[1] + resi_ratio4 * kFactors[2]) * d4start +
            (kFactors[0] + ratio2 * kFactors[1] + ratio4 * kFactors[2]) * d4end) *
               scale * scale * scale;
}

template <std::floating_point T = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto quinerpCoeffs(T ratio, T scale = 1.0) noexcept -> std::array<T, 6> {
    constexpr std::array<T, 5> kFactors{
        -(1.0 / 6.0), (1.0 / 6.0), 7.0 / 360.0, -(1.0 / 36.0), 1.0 / 120.0
    };
    const T ratio3{ratio * ratio * ratio};
    const T ratio5{ratio3 * ratio * ratio};
    const T resi_ratio{1.0 - ratio};
    const T resi_ratio3{resi_ratio * resi_ratio * resi_ratio};
    const T resi_ratio5{resi_ratio3 * resi_ratio * resi_ratio};
    const T scale2{scale * scale};
    const T scale4{scale2 * scale2};
    return std::array<T, 6>{
        resi_ratio,
        ratio,
        (resi_ratio * kFactors[0] + resi_ratio3 * kFactors[1]) * scale2,
        (ratio * kFactors[0] + ratio3 * kFactors[1]) * scale2,
        (resi_ratio * kFactors[2] + resi_ratio3 * kFactors[3] + resi_ratio5 * kFactors[4]) * scale4,
        (ratio * kFactors[2] + ratio3 * kFactors[3] + ratio5 * kFactors[4]) * scale4
    };
}

template <std::floating_point T = double>
[[using gnu: const, always_inline, hot]] [[nodiscard]]
inline constexpr auto quinerpdCoeffs(T ratio, T scale = 1.0) noexcept -> std::array<T, 6> {
    constexpr std::array<T, 5> kFactors{-(1.0 / 6.0), 0.5, 7.0 / 360.0, -(1.0 / 12.0), 1.0 / 24.0};
    const T ratio2{ratio * ratio};
    const T ratio4{ratio2 * ratio2};
    const T resi_ratio{1.0 - ratio};
    const T resi_ratio2{resi_ratio * resi_ratio};
    const T resi_ratio4{resi_ratio2 * resi_ratio2};
    const T reci_scale{1.0 / scale};
    const T scale3{scale * scale * scale};
    return std::array<T, 6>{
        -reci_scale,
        reci_scale,
        -(kFactors[0] + resi_ratio2 * kFactors[1]) * scale,
        (kFactors[0] + ratio2 * kFactors[1]) * scale,
        -(kFactors[2] + resi_ratio2 * kFactors[3] + resi_ratio4 * kFactors[4]) * scale3,
        (kFactors[2] + ratio2 * kFactors[3] + ratio4 * kFactors[4]) * scale3
    };
}

} // namespace boyle::math
