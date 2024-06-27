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
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/macros.hpp"
#include "boyle/kinetics/dualism.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class [[nodiscard]] Border2 {
    friend class boost::serialization::access;

  public:
    explicit Border2(
        std::uint64_t c_id, ::boyle::kinetics::Chirality c_chirality,
        std::vector<::boyle::math::Vec2<T>> c_bound_points
    ) noexcept
        : id{c_id}, chirality{c_chirality}, bound_points{std::move(c_bound_points)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(Border2);
    virtual ~Border2() noexcept = 0;

    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Chirality chirality{::boyle::kinetics::Chirality::LEFT};
    std::vector<::boyle::math::Vec2<T>> bound_points{};

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & id;
        archive & chirality;
        archive & bound_points;
        return;
    }
};

template <std::floating_point T>
inline Border2<T>::~Border2() noexcept = default;

template <std::floating_point T>
class [[nodiscard]] HardBorder2 final : public Border2<T> {
    friend class boost::serialization::access;

  public:
    explicit HardBorder2(
        std::uint64_t c_id, ::boyle::kinetics::Chirality c_chirality,
        std::vector<::boyle::math::Vec2<T>> c_bound_points
    ) noexcept
        : Border2<T>{c_id, c_chirality, std::move(c_bound_points)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(HardBorder2);
    ~HardBorder2() noexcept override = default;

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive& boost::serialization::base_object<Border2<T>>(*this);
        return;
    }
};

template <std::floating_point T>
class [[nodiscard]] SoftBorder2 final : public Border2<T> {
    friend class boost::serialization::access;

  public:
    explicit SoftBorder2(
        std::uint64_t c_id, ::boyle::kinetics::Chirality c_chirality,
        std::vector<::boyle::math::Vec2<T>> c_bound_points, T c_linear_weight, T c_quadratic_weight
    ) noexcept
        : Border2<T>{c_id, c_chirality, std::move(c_bound_points)}, linear_weight{c_linear_weight},
          quadratic_weight{c_quadratic_weight} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(SoftBorder2);
    ~SoftBorder2() noexcept override = default;

    T linear_weight{0.0};
    T quadratic_weight{0.0};

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive& boost::serialization::base_object<Border2<T>>(*this);
        archive & linear_weight;
        archive & quadratic_weight;
        return;
    }
};

using HardBorder2f = HardBorder2<float>;
using HardBorder2d = HardBorder2<double>;

using SoftBorder2f = SoftBorder2<float>;
using SoftBorder2d = SoftBorder2<double>;

} // namespace boyle::kinetics
