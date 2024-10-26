/**
 * @file piecewise_cubic_curve.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-01
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "fmt/format.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/cubic_interpolation.hpp"
#include "boyle/math/duplet.hpp"
#include "boyle/math/functions/piecewise_cubic_function1.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"
#include "boyle/math/vec3.hpp"

namespace boyle::math {

template <VecArithmetic T, std::floating_point U = typename T::value_type>
class [[nodiscard]] PiecewiseCubicCurve final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;
    using BoundaryMode = typename PiecewiseCubicFunction1<value_type, param_type>::BoundaryMode;

    static constexpr param_type kDuplicateCriterion{1E-8};

    PiecewiseCubicCurve() noexcept = default;
    PiecewiseCubicCurve(const PiecewiseCubicCurve& other) noexcept = default;
    auto operator=(const PiecewiseCubicCurve& other) noexcept -> PiecewiseCubicCurve& = default;
    PiecewiseCubicCurve(PiecewiseCubicCurve&& other) noexcept = default;
    auto operator=(PiecewiseCubicCurve&& other) noexcept -> PiecewiseCubicCurve& = default;
    ~PiecewiseCubicCurve() noexcept = default;

    [[using gnu: always_inline]]
    explicit PiecewiseCubicCurve(std::vector<value_type> anchor_points, param_type s0 = 0.0)
        : PiecewiseCubicCurve(
              std::move(anchor_points), BoundaryMode{2, value_type{0.0, 0.0}},
              BoundaryMode{2, value_type{0.0, 0.0}}, s0
          ) {}

    [[using gnu: ]]
    explicit PiecewiseCubicCurve(
        std::vector<value_type> anchor_points, BoundaryMode b0, BoundaryMode bf, param_type s0 = 0.0
    ) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
        if (anchor_points.size() < 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! sizes of anchor_points must be greater than 2: "
                "anchor_points.size() = {0:d}",
                anchor_points.size()
            ));
        }
#endif
        const std::size_t size{anchor_points.size()};
        std::vector<param_type> temp_arc_lengths(size);
        temp_arc_lengths[0] = s0;
        for (std::size_t i{1}; i < size; ++i) {
            temp_arc_lengths[i] =
                temp_arc_lengths[i - 1] + anchor_points[i].euclideanTo(anchor_points[i - 1]);
        }
        const PiecewiseCubicFunction1<value_type, param_type> temp_vec_of_s{
            temp_arc_lengths, anchor_points, b0, bf
        };

        const std::vector<value_type>& ddys = temp_vec_of_s.ddys();
        std::vector<param_type> arc_lengths(size);
        arc_lengths[0] = s0;
        for (std::size_t i{1}; i < size; ++i) {
            arc_lengths[i] =
                arc_lengths[i - 1] + ::boyle::math::calcArcLength(
                                         anchor_points[i - 1], anchor_points[i], ddys[i - 1],
                                         ddys[i], temp_arc_lengths[i] - temp_arc_lengths[i - 1]
                                     );
        }
        m_vec_of_s = PiecewiseCubicFunction1<value_type, param_type>{
            std::move(arc_lengths), std::move(anchor_points), b0, bf
        };
    }

    [[using gnu: ]]
    explicit PiecewiseCubicCurve(
        [[maybe_unused]] periodic_tag tag, std::vector<value_type> anchor_points,
        param_type s0 = 0.0
    ) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
        if (anchor_points.size() < 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! sizes of anchor_points must be greater than 2: "
                "anchor_points.size() = {0:d}",
                anchor_points.size()
            ));
        }
#endif
        const std::size_t size{anchor_points.size()};
        std::vector<param_type> temp_arc_lengths(size);
        temp_arc_lengths[0] = s0;
        for (std::size_t i{1}; i < size; ++i) {
            temp_arc_lengths[i] =
                temp_arc_lengths[i - 1] + anchor_points[i].euclideanTo(anchor_points[i - 1]);
        }
        const PiecewiseCubicFunction1<value_type, param_type> temp_vec_of_s{
            tag, temp_arc_lengths, anchor_points
        };

        const std::vector<value_type>& ddys = temp_vec_of_s.ddys();
        std::vector<param_type> arc_lengths(size);
        arc_lengths[0] = s0;
        for (std::size_t i{1}; i < size; ++i) {
            arc_lengths[i] =
                arc_lengths[i - 1] + ::boyle::math::calcArcLength(
                                         anchor_points[i - 1], anchor_points[i], ddys[i - 1],
                                         ddys[i], temp_arc_lengths[i] - temp_arc_lengths[i - 1]
                                     );
        }
        m_vec_of_s = PiecewiseCubicFunction1<value_type, param_type>{
            tag, std::move(arc_lengths), std::move(anchor_points)
        };
    }

    [[using gnu: pure, always_inline]]
    auto eval(param_type s) const noexcept -> value_type {
        return m_vec_of_s.eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto eval(param_type s, param_type l) const noexcept -> value_type
        requires InstanceOfTemplate<value_type, Vec2>
    {
        return eval(s) + l * orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    auto eval(SlDuplet<param_type> sl) const noexcept -> value_type
        requires InstanceOfTemplate<value_type, Vec2>
    {
        return eval(sl.s) + sl.l * orthonormal(sl.s);
    }

    [[using gnu: pure, always_inline]]
    auto tangent(param_type s) const noexcept -> value_type {
        return m_vec_of_s.derivative(s).normalized();
    }

    [[using gnu: pure, always_inline]]
    auto orthonormal(param_type s) const noexcept -> value_type
        requires InstanceOfTemplate<value_type, Vec2>
    {
        return tangent(s).selfRotateHalfPi();
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto curvature(param_type s) const noexcept -> param_type {
        const std::vector<param_type>& arc_lengths{arcLengths()};
        const std::vector<value_type>& anchor_points{anchorPoints()};
        const std::vector<value_type>& ddys{m_vec_of_s.ddys()};
        const std::size_t size{arc_lengths.size()};
        const std::size_t pos =
            nearestUpperElement(
                std::ranges::subrange{arc_lengths.cbegin(), arc_lengths.cend()}, s
            ) -
            arc_lengths.cbegin();
        if (pos == 0) {
            return 0.0;
        }
        if (pos == size) {
            return 0.0;
        }
        const param_type h{arc_lengths[pos] - arc_lengths[pos - 1]};
        const param_type ratio{(s - arc_lengths[pos - 1]) / h};
        const value_type derivative{
            cuberpd(anchor_points[pos - 1], anchor_points[pos], ddys[pos - 1], ddys[pos], ratio, h)
        };
        const value_type derivative2{lerp(ddys[pos - 1], ddys[pos], ratio)};
        const param_type derivative_norm{derivative.norm()};
        return derivative.crossNorm(derivative2) /
               (derivative_norm * derivative_norm * derivative_norm);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(value_type point) const noexcept -> SlDuplet<param_type>
        requires InstanceOfTemplate<value_type, Vec2>
    {
        const std::vector<param_type>& arc_lengths{arcLengths()};
        const std::vector<value_type>& anchor_points{anchorPoints()};
        const std::vector<value_type>& ddys{m_vec_of_s.ddys()};
        const std::size_t pos =
            nearestUpperElement(
                std::ranges::subrange{anchor_points.cbegin(), anchor_points.cend()}, point
            ) -
            anchor_points.cbegin();
        value_type derivative, derivative2, r;
        param_type h, derivative_norm, s, l;
        if (pos == 0) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = cuberpd(anchor_points[0], anchor_points[1], ddys[0], ddys[1], 0.0, h);
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.crossNorm(r) / derivative_norm;
        } else if (pos == anchor_points.size()) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = cuberpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.crossNorm(r) / derivative_norm;
        } else {
            const value_type& lower_anchor_point{anchor_points[pos - 1]};
            const value_type& upper_anchor_point{anchor_points[pos]};
            const value_type& lower_derivative2{ddys[pos - 1]};
            const value_type& upper_derivative2{ddys[pos]};
            param_type ratio = 0.5;
            h = m_vec_of_s.ts()[pos] - m_vec_of_s.ts()[pos - 1];
            r = point - cuberp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, ratio, h
                        );
            derivative = cuberpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2, ratio,
                h
            );
            derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.normSqr()) * h);
                r = point - cuberp(
                                lower_anchor_point, upper_anchor_point, lower_derivative2,
                                upper_derivative2, ratio, h
                            );
                derivative = cuberpd(
                    lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                    ratio, h
                );
                derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            }
            s = lerp(arc_lengths[pos - 1], arc_lengths[pos], ratio);
            l = derivative.crossNorm(r) / derivative.norm();
        }
        return SlDuplet<param_type>{s, l};
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(value_type point, param_type start_s, param_type end_s) const noexcept
        -> SlDuplet<param_type>
        requires InstanceOfTemplate<value_type, Vec2>
    {
        if (start_s > end_s) {
            std::swap(start_s, end_s);
        }
        const std::vector<param_type>& arc_lengths{arcLengths()};
        const std::vector<value_type>& anchor_points{anchorPoints()};
        const std::vector<value_type>& ddys{m_vec_of_s.ddys()};
        const std::size_t istart =
            nearestUpperElement(
                std::ranges::subrange{arc_lengths.cbegin(), arc_lengths.cend()}, start_s
            ) -
            arc_lengths.cbegin() - 1;
        const std::size_t iend =
            nearestUpperElement(
                std::ranges::subrange{arc_lengths.cbegin(), arc_lengths.cend()}, end_s
            ) -
            arc_lengths.cbegin();
        const std::size_t pos =
            nearestUpperElement(
                std::ranges::subrange{
                    anchor_points.cbegin() + istart, anchor_points.cbegin() + iend
                },
                point
            ) -
            anchor_points.cbegin();
        value_type derivative, derivative2, r;
        param_type h, derivative_norm, s, l;
        if (pos == istart) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = cuberpd(anchor_points[0], anchor_points[1], ddys[0], ddys[1], 0.0, h);
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.crossNorm(r) / derivative_norm;
        } else if (pos == iend) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = cuberpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.crossNorm(r) / derivative_norm;
        } else {
            const value_type& lower_anchor_point{anchor_points[pos - 1]};
            const value_type& upper_anchor_point{anchor_points[pos]};
            const value_type& lower_derivative2{ddys[pos - 1]};
            const value_type& upper_derivative2{ddys[pos]};
            param_type ratio = 0.5;
            h = m_vec_of_s.ts()[pos] - m_vec_of_s.ts()[pos - 1];
            r = point - cuberp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, ratio, h
                        );
            derivative = cuberpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2, ratio,
                h
            );
            derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.normSqr()) * h);
                r = point - cuberp(
                                lower_anchor_point, upper_anchor_point, lower_derivative2,
                                upper_derivative2, ratio, h
                            );
                derivative = cuberpd(
                    lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                    ratio, h
                );
                derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            }
            s = lerp(arc_lengths[pos - 1], arc_lengths[pos], ratio);
            l = derivative.crossNorm(r) / derivative.norm();
        }
        return SlDuplet<param_type>{s, l};
    }

    [[using gnu: pure, always_inline]]
    auto operator()(param_type s) const noexcept -> value_type {
        return eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(param_type s, param_type l) const noexcept -> value_type
        requires InstanceOfTemplate<value_type, Vec2>
    {
        return eval(s, l);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(SlDuplet<param_type> sl) const noexcept -> value_type
        requires InstanceOfTemplate<value_type, Vec2>
    {
        return eval(sl);
    }

    [[using gnu: pure, always_inline]]
    auto minS() const noexcept -> param_type {
        return m_vec_of_s.minT();
    }

    [[using gnu: pure, always_inline]]
    auto maxS() const noexcept -> param_type {
        return m_vec_of_s.maxT();
    }

    [[using gnu: pure, always_inline]]
    auto front() const noexcept -> value_type {
        return m_vec_of_s.ys().front();
    }

    [[using gnu: pure, always_inline]]
    auto back() const noexcept -> value_type {
        return m_vec_of_s.ys().back();
    }

    [[using gnu: pure, always_inline]]
    auto arcLengths() const noexcept -> const std::vector<param_type>& {
        return m_vec_of_s.ts();
    }

    [[using gnu: pure, always_inline]]
    auto anchorPoints() const noexcept -> const std::vector<value_type>& {
        return m_vec_of_s.ys();
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_vec_of_s;
        return;
    }

    PiecewiseCubicFunction1<value_type, param_type> m_vec_of_s{};
};

using PiecewiseCubicCurve2f = PiecewiseCubicCurve<Vec2f>;
using PiecewiseCubicCurve2d = PiecewiseCubicCurve<Vec2d>;

using PiecewiseCubicCurve3f = PiecewiseCubicCurve<Vec3f>;
using PiecewiseCubicCurve3d = PiecewiseCubicCurve<Vec3d>;

} // namespace boyle::math
