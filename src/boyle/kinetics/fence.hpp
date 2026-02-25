/**
 * @file fence.hpp
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

namespace boyle::kinetics {

template <std::floating_point T>
struct HardFence final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Actio actio{::boyle::kinetics::Actio::BLOCKING};
    std::pmr::vector<value_type> bound_ts;
    std::pmr::vector<value_type> bound_ss;
};

template <std::floating_point T>
struct SoftFence final {
    using value_type = T;
    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Actio actio{::boyle::kinetics::Actio::BLOCKING};
    std::pmr::vector<value_type> bound_ts;
    std::pmr::vector<value_type> bound_ss;
    T linear_weight{0.0};
    T quadratic_weight{0.0};
};

using HardFence1s = HardFence<float>;

using HardFence1d = HardFence<double>;

using SoftFence1s = SoftFence<float>;

using SoftFence1d = SoftFence<double>;

} // namespace boyle::kinetics

namespace boost::serialization {

template <std::floating_point T>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::kinetics::HardFence<T>& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.id;
    archive & obj.actio;
    archive & obj.bound_ts;
    archive & obj.bound_ss;
    return;
}

template <std::floating_point T>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::kinetics::SoftFence<T>& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.id;
    archive & obj.actio;
    archive & obj.bound_ts;
    archive & obj.bound_ss;
    archive & obj.linear_weight;
    archive & obj.quadratic_weight;
    return;
}

} // namespace boost::serialization
