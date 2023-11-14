/**
 * @file border1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "common/dualism.hpp"
#include "math/functions/piecewise_functions/piecewise_linear_function1.hpp"

namespace tiny_pnc {
namespace kinetics {

template <tiny_pnc::common::Stiffness SF, typename T>
struct Border1;

template <typename T>
struct [[nodiscard]] Border1<tiny_pnc::common::Stiffness::HARD, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    std::uint64_t id;
    tiny_pnc::common::Chirality chirality;
    math::PiecewiseLinearFunction1<T> bound_line;
};

template <typename T>
struct [[nodiscard]] Border1<tiny_pnc::common::Stiffness::SOFT, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    std::uint64_t id;
    tiny_pnc::common::Chirality chirality;
    math::PiecewiseLinearFunction1<T> bound_line;
    T linear_weight;
    T quadratic_weight;
};

template <typename T>
using HardBorder1 = Border1<tiny_pnc::common::Stiffness::HARD, T>;

template <typename T>
using SoftBorder1 = Border1<tiny_pnc::common::Stiffness::SOFT, T>;

using HardBorder1f = HardBorder1<float>;
using HardBorder1d = HardBorder1<double>;

using SoftBorder1f = SoftBorder1<float>;
using SoftBorder1d = SoftBorder1<double>;

} // namespace kinetics
} // namespace tiny_pnc

namespace boost {
namespace serialization {

template <typename Archive, typename T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::kinetics::Border1<tiny_pnc::common::Stiffness::HARD, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.chirality;
    ar& obj.bound_line;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::kinetics::Border1<tiny_pnc::common::Stiffness::SOFT, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.chirality;
    ar& obj.bound_line;
    ar& obj.linear_weight;
    ar& obj.quadratic_weight;
    return;
}

} // namespace serialization
} // namespace boost
