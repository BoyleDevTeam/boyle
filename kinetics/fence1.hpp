/**
 * @file fence1.hpp
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
struct Fence1;

template <typename T>
struct [[nodiscard]] Fence1<tiny_pnc::common::Stiffness::HARD, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    std::uint64_t id;
    tiny_pnc::common::Actio actio;
    math::PiecewiseLinearFunction1<T> bound_line;
};

template <typename T>
struct [[nodiscard]] Fence1<tiny_pnc::common::Stiffness::SOFT, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    std::uint64_t id;
    tiny_pnc::common::Actio actio;
    math::PiecewiseLinearFunction1<T> bound_line;
    T linear_weight;
    T quadratic_weight;
};

template <typename T>
using HardFence1 = Fence1<tiny_pnc::common::Stiffness::HARD, T>;

template <typename T>
using SoftFence1 = Fence1<tiny_pnc::common::Stiffness::SOFT, T>;

using HardFence1f = HardFence1<float>;
using HardFence1d = HardFence1<double>;

using SoftFence1f = SoftFence1<float>;
using SoftFence1d = SoftFence1<double>;

} // namespace kinetics
} // namespace tiny_pnc

namespace boost {
namespace serialization {

template <typename Archive, typename T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::kinetics::Fence1<tiny_pnc::common::Stiffness::HARD, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.actio;
    ar& obj.bound_line;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::kinetics::Fence1<tiny_pnc::common::Stiffness::SOFT, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.actio;
    ar& obj.bound_line;
    ar& obj.linear_weight;
    ar& obj.quadratic_weight;
    return;
}

} // namespace serialization
} // namespace boost
