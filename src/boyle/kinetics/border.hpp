/**
 * @file border.hpp
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

#include <cstdint>
#include <limits>
#include <vector>

#include "boost/serialization/vector.hpp"

#include "boyle/kinetics/dualism.hpp"
#include "boyle/math/dense/vec2.hpp"

namespace boyle::kinetics {

template <::boyle::math::Vec2Arithmetic T>
struct HardBorder final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Chirality chirality{::boyle::kinetics::Chirality::LEFT};
    std::pmr::vector<value_type> bound_points;
};

template <::boyle::math::Vec2Arithmetic T>
struct SoftBorder final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Chirality chirality{::boyle::kinetics::Chirality::LEFT};
    std::pmr::vector<value_type> bound_points;
    typename value_type::value_type linear_weight{0.0};
    typename value_type::value_type quadratic_weight{0.0};
};

using HardBorder2s = HardBorder<::boyle::math::Vec2s>;

using HardBorder2d = HardBorder<::boyle::math::Vec2d>;

using SoftBorder2s = SoftBorder<::boyle::math::Vec2s>;

using SoftBorder2d = SoftBorder<::boyle::math::Vec2d>;

} // namespace boyle::kinetics

namespace boost::serialization {

template <boyle::math::Vec2Arithmetic T>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::kinetics::HardBorder<T>& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.id;
    archive & obj.chirality;
    archive & obj.bound_points;
    return;
}

template <boyle::math::Vec2Arithmetic T>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::kinetics::SoftBorder<T>& obj,
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
