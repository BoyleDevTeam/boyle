/**
 * @file piecewise_linear_curve.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-28
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "fmt/format.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/duplet.hpp"
#include "boyle/math/functions/piecewise_linear_function1.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"
#include "boyle/math/vec3.hpp"

namespace boyle::math {

template <VecArithmetic T, std::floating_point U = typename T::value_type>
class [[nodiscard]] PiecewiseLinearCurve final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr value_type kDuplicateCriterion{1E-8};

    PiecewiseLinearCurve() noexcept = default;
    PiecewiseLinearCurve(const PiecewiseLinearCurve& other) noexcept = default;
    auto operator=(const PiecewiseLinearCurve& other) noexcept -> PiecewiseLinearCurve& = default;
    PiecewiseLinearCurve(PiecewiseLinearCurve&& other) noexcept = default;
    auto operator=(PiecewiseLinearCurve&& other) noexcept -> PiecewiseLinearCurve& = default;
    ~PiecewiseLinearCurve() noexcept = default;

    [[using gnu: ]]
    explicit PiecewiseLinearCurve(
        std::vector<value_type> anchor_points, param_type s0 = 0.0
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
        const std::size_t size = anchor_points.size();
        std::vector<param_type> arc_lengths(size);
        arc_lengths[0] = s0;
        for (std::size_t i{1}; i < size; ++i) {
            arc_lengths[i] =
                arc_lengths[i - 1] + anchor_points[i].euclideanTo(anchor_points[i - 1]);
        }
        m_vec_of_s = PiecewiseLinearFunction1<value_type, param_type>{
            std::move(arc_lengths), std::move(anchor_points)
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
    auto inverse(value_type point) const noexcept -> SlDuplet<param_type>
        requires InstanceOfTemplate<value_type, Vec2>
    {
        const std::vector<param_type>& arc_lengths = arcLengths();
        const std::vector<value_type>& anchor_points = anchorPoints();
        const std::size_t pos =
            nearestUpperElement(
                std::ranges::subrange{anchor_points.cbegin(), anchor_points.cend()}, point
            ) -
            anchor_points.cbegin();
        value_type diff, r;
        param_type diff_norm, s, l;
        if (pos == 0) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_norm = diff.norm();
            s = arc_lengths.front() + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
        } else if (pos == anchor_points.size()) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_norm = diff.norm();
            s = arc_lengths.back() + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_norm = diff.norm();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
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
        const std::vector<param_type>& arc_lengths = arcLengths();
        const std::vector<value_type>& anchor_points = anchorPoints();
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
        value_type diff, r;
        param_type diff_norm, s, l;
        if (pos == istart) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_norm = diff.norm();
            s = arc_lengths.front() + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
        } else if (pos == iend) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_norm = diff.norm();
            s = arc_lengths.back() + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_norm = diff.norm();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_norm;
            l = diff.crossNorm(r) / diff_norm;
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

    PiecewiseLinearFunction1<value_type, param_type> m_vec_of_s{};
};

using PiecewiseLinearCurve2f = PiecewiseLinearCurve<Vec2f>;
using PiecewiseLinearCurve2d = PiecewiseLinearCurve<Vec2d>;

using PiecewiseLinearCurve3f = PiecewiseLinearCurve<Vec3f>;
using PiecewiseLinearCurve3d = PiecewiseLinearCurve<Vec3d>;

} // namespace boyle::math
