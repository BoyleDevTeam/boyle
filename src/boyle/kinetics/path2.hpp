/**
 * @file path2.hpp
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

#include <array>
#include <concepts>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "boyle/math/curves/piecewise_quintic_curve.hpp"
#include "boyle/math/dense/vec2.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class [[nodiscard]] Path2 final {
    friend class boost::serialization::access;

  public:
    using BoundaryMode =
        typename ::boyle::math::PiecewiseQuinticCurve<::boyle::math::Vec2<T>>::BoundaryMode;

    Path2() noexcept = default;
    Path2(const Path2& other) noexcept = default;
    auto operator=(const Path2& other) noexcept -> Path2& = default;
    Path2(Path2&& other) noexcept = default;
    auto operator=(Path2&& other) noexcept -> Path2& = default;
    ~Path2() noexcept = default;

    [[using gnu: always_inline]]
    explicit Path2(std::vector<::boyle::math::Vec2<T>> anchor_points, T s0 = 0.0)
        : Path2{
              std::move(anchor_points),
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, ::boyle::math::Vec2<T>{0.0, 0.0}},
                  BoundaryMode{4, ::boyle::math::Vec2<T>{0.0, 0.0}}
              },
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, ::boyle::math::Vec2<T>{0.0, 0.0}},
                  BoundaryMode{4, ::boyle::math::Vec2<T>{0.0, 0.0}}
              },
              s0
          } {}

    [[using gnu: always_inline]]
    explicit Path2(
        std::vector<::boyle::math::Vec2<T>> anchor_points, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf, T s0 = 0.0
    )
        : m_curve{anchor_points, b0, bf, s0} {}

    [[using gnu: pure, always_inline]]
    auto operator()(T s) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve(s);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(T s, T l) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve(s, l);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(::boyle::math::SlDuplet<T> sl) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve(sl);
    }

    [[using gnu: pure, always_inline]]
    auto inverse(::boyle::math::Vec2<T> point) const noexcept -> ::boyle::math::SlDuplet<T> {
        return m_curve.inverse(point);
    }

    [[using gnu: pure, always_inline]]
    auto inverse(::boyle::math::Vec2<T> point, T start_s, T end_s) const noexcept
        -> ::boyle::math::SlDuplet<T> {
        return m_curve.inverse(point, start_s, end_s);
    }

    [[using gnu: pure, always_inline]]
    auto tangent(T s) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve.tangent(s);
    }

    [[using gnu: pure, always_inline]]
    auto normal(T s) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve.normal(s);
    }

    [[using gnu: pure, always_inline]]
    auto curvature(T s) const noexcept -> T {
        return m_curve.curvature(s);
    }

    [[using gnu: pure, always_inline]]
    auto minS() const noexcept -> T {
        return m_curve.minS();
    }

    [[using gnu: pure, always_inline]]
    auto maxS() const noexcept -> T {
        return m_curve.maxS();
    }

    [[using gnu: pure, always_inline]]
    auto arcLengths() const noexcept -> std::span<const T> {
        return m_curve.arcLengths();
    }

    [[using gnu: pure, always_inline]]
    auto anchorPoints() const noexcept -> std::span<const ::boyle::math::Vec2<T>> {
        return m_curve.anchorPoints();
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive) noexcept -> void {
        archive & m_curve;
        return;
    }

    ::boyle::math::PiecewiseQuinticCurve<::boyle::math::Vec2<T>> m_curve{};
};

using Path2f = Path2<float>;
using Path2d = Path2<double>;

} // namespace boyle::kinetics
