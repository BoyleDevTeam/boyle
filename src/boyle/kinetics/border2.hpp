/**
 * @file border2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <cstdint>
#include <limits>
#include <vector>

#include "boost/serialization/vector.hpp"

#include "boyle/kinetics/dualism.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
struct [[nodiscard]] HardBorder2 final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Chirality chirality{::boyle::kinetics::Chirality::LEFT};
    std::vector<::boyle::math::Vec2<value_type>> bound_points{};
};

template <std::floating_point T>
struct [[nodiscard]] SoftBorder2 final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Chirality chirality{::boyle::kinetics::Chirality::LEFT};
    std::vector<::boyle::math::Vec2<value_type>> bound_points{};
    value_type linear_weight{0.0};
    value_type quadratic_weight{0.0};
};

using HardBorder2f = HardBorder2<float>;
using HardBorder2d = HardBorder2<double>;

using SoftBorder2f = SoftBorder2<float>;
using SoftBorder2d = SoftBorder2<double>;

} // namespace boyle::kinetics

namespace boost::serialization {

[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, boyle::math::InstanceOfTemplate<boyle::kinetics::HardBorder2> auto& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.id;
    archive & obj.chirality;
    archive & obj.bound_points;
    return;
}

[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, boyle::math::InstanceOfTemplate<boyle::kinetics::SoftBorder2> auto& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.id;
    archive & obj.chirality;
    archive & obj.bound_points;
    archive & obj.linear_weight;
    archive & obj.quadratic_weight;
    return;
}

} // namespace boost::serialization
