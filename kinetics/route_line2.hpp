/**
 * @file route_line2.hpp
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
#include <type_traits>
#include <vector>

#include "boost/serialization/access.hpp"

#include "common/utils/macros.hpp"
#include "math/curves/piecewise_curves/piecewise_quintic_curve2.hpp"

namespace tiny_pnc {
namespace kinetics {

template <typename T>
class [[nodiscard]] RouteLine2 final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using BoundaryMode = typename tiny_pnc::math::PiecewiseQuinticCurve2<T>::BoundaryMode;

    [[using gnu: always_inline]] explicit RouteLine2(
        std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>> anchor_points, T s0 = 0.0
    )
        : RouteLine2(
              std::move(anchor_points),
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>{0.0, 0.0}},
                  BoundaryMode{4, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>{0.0, 0.0}}},
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>{0.0, 0.0}},
                  BoundaryMode{4, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>{0.0, 0.0}}},
              s0
          ) {}

    [[using gnu: always_inline]] explicit RouteLine2(
        std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>> anchor_points,
        std::array<BoundaryMode, 2> b0, std::array<BoundaryMode, 2> bf, T s0 = 0.0
    )
        : curve_(anchor_points, b0, bf, s0) {}

    ENABLE_IMPLICIT_CONSTRUCTORS(RouteLine2);

    ~RouteLine2() noexcept = default;

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>
    operator()(T s) const noexcept {
        return curve_(s);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>
    operator()(T s, T l) const noexcept {
        return curve_(s, l);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>
    operator()(tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T> sl) const noexcept {
        return curve_(sl);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T> inverse(
        tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> point
    ) const noexcept {
        return curve_.inverse(point);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T> inverse(
        tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> point, T start_s, T end_s
    ) const noexcept {
        return curve_.inverse(point, start_s, end_s);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> tangent(T s) const noexcept {
        return curve_.tangent(s);
    }

    [[using gnu: pure, always_inline]]
    tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> orthnormal(T s) const noexcept {
        return curve_.orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    T curvature(T s) {
        return curve_.curvature(s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>>
    operator()(const std::vector<T>& s) const noexcept {
        return curve_(s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>> tangent(
        const std::vector<T>& s
    ) const noexcept {
        return curve_.tangent(s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>> orthonormal(
        const std::vector<T>& s
    ) const noexcept {
        return curve_.orthonormal(s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>>
    operator()(const std::vector<tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T>> sls
    ) const noexcept {
        return curve_(sls);
    }

  private:
    [[using gnu: pure,
      always_inline]] explicit RouteLine2(tiny_pnc::math::PiecewiseQuinticCurve2<T> curve) noexcept
        : curve_(std::move(curve)) {}

    tiny_pnc::math::PiecewiseQuinticCurve2<T> curve_;
};

using RouteLine2f = RouteLine2<float>;
using RouteLine2d = RouteLine2<double>;

} // namespace kinetics
} // namespace tiny_pnc
