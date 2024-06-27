/**
 * @file route_line2.hpp
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

#include "boyle/common/utils/macros.hpp"
#include "boyle/math/curves/piecewise_curves/piecewise_quintic_curve2.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class [[nodiscard]] RouteLine2 final {
    friend class boost::serialization::access;

  public:
    using BoundaryMode = typename ::boyle::math::PiecewiseQuinticCurve2<T>::BoundaryMode;

    [[using gnu: always_inline]]
    explicit RouteLine2(std::vector<::boyle::math::Vec2<T>> anchor_points, T s0 = 0.0)
        : RouteLine2{
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
    explicit RouteLine2(
        std::vector<::boyle::math::Vec2<T>> anchor_points, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf, T s0 = 0.0
    )
        : m_curve{anchor_points, b0, bf, s0} {}

    ENABLE_IMPLICIT_CONSTRUCTORS(RouteLine2);
    ~RouteLine2() noexcept = default;

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
    auto orthnormal(T s) const noexcept -> ::boyle::math::Vec2<T> {
        return m_curve.orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    auto curvature(T s) -> T {
        return m_curve.curvature(s);
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto operator()(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<::boyle::math::Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        return m_curve(first, last);
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto tangent(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<::boyle::math::Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        return m_curve.tangent(first, last);
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto orthonormal(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<::boyle::math::Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        return m_curve.orthonormal(first, last);
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto operator()(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<::boyle::math::Vec2<T>>
        requires std::same_as<
            ::boyle::math::SlDuplet<T>, typename std::iterator_traits<ForwardIt>::value_type>
    {
        return m_curve(first, last);
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
    auto arcLengths() const noexcept -> const std::vector<T>& {
        return m_curve.arcLengths();
    }

    [[using gnu: pure, always_inline]]
    auto anchorPoints() const noexcept -> const std::vector<::boyle::math::Vec2<T>>& {
        return m_curve.anchorPoints();
    }

  private:
    [[using gnu: always_inline]]
    explicit RouteLine2(::boyle::math::PiecewiseQuinticCurve2<T> curve) noexcept
        : m_curve{std::move(curve)} {}

    [[using gnu: always_inline]]
    auto serialize(auto& archive) noexcept -> void {
        archive & m_curve;
        return;
    }

    ::boyle::math::PiecewiseQuinticCurve2<T> m_curve{};
};

using RouteLine2f = RouteLine2<float>;
using RouteLine2d = RouteLine2<double>;

} // namespace boyle::kinetics
