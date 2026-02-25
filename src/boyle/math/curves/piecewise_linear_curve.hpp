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
#include <format>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "boyle/common/utils/aligned_allocator.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/curves/sl.hpp"
#include "boyle/math/curves/slv.hpp"
#include "boyle/math/dense/vec2.hpp"
#include "boyle/math/dense/vec3.hpp"
#include "boyle/math/functions/piecewise_linear_function.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <VecArithmetic T, Allocatory Alloc = ::boyle::common::AlignedAllocator<T, alignof(T)>>
    requires std::floating_point<typename T::value_type>
class PiecewiseLinearCurve final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = typename T::value_type;
    using size_type = std::size_t;
    using allocator_type = Alloc;

    static constexpr value_type kDuplicateCriterion{1E-8};

    PiecewiseLinearCurve() noexcept = default;
    PiecewiseLinearCurve(const PiecewiseLinearCurve& other) = default;
    auto operator=(const PiecewiseLinearCurve& other) -> PiecewiseLinearCurve& = default;
    PiecewiseLinearCurve(PiecewiseLinearCurve&& other) noexcept = default;
    auto operator=(PiecewiseLinearCurve&& other) noexcept -> PiecewiseLinearCurve& = default;
    ~PiecewiseLinearCurve() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_vec_of_s.get_allocator();
    }

    template <std::ranges::input_range R = std::initializer_list<value_type>>
    [[using gnu: ]]
    explicit PiecewiseLinearCurve(
        R&& anchor_points, param_type s0 = 0.0, const allocator_type& alloc = {}
    )
        requires std::same_as<std::ranges::range_value_t<R>, value_type>
    {
#if BOYLE_CHECK_PARAMS == 1
        if (anchor_points.size() < 2) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! sizes of anchor_points must be greater than 2: "
                    "anchor_points.size() = {0:d}",
                    anchor_points.size()
                )
            );
        }
#endif
        const size_type size{anchor_points.size()};
        std::vector<param_type, param_allocator_type> arc_lengths(size, alloc);
        arc_lengths[0] = s0;
        for (size_type i{1}; i < size; ++i) {
            arc_lengths[i] =
                arc_lengths[i - 1] + anchor_points[i].euclideanTo(anchor_points[i - 1]);
        }
        m_vec_of_s =
            PiecewiseLinearFunction<value_type, allocator_type>(arc_lengths, anchor_points, alloc);
    }

    [[using gnu: pure, always_inline]]
    auto eval(param_type s) const noexcept -> value_type {
        return m_vec_of_s.eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto eval(param_type s, param_type l) const noexcept -> value_type
        requires Vec2Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{arc_lengths.size()};
        const size_type pos = nearestUpperElement(arc_lengths, s) - arc_lengths.begin();
        if (pos < 2) {
            const param_type ratio{(s - arc_lengths[0]) / (arc_lengths[1] - arc_lengths[0])};
            const value_type val{lerp(anchor_points[0], anchor_points[1], ratio)};
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return val + normal * l;
        }
        if (pos == size) {
            const param_type ratio{
                (s - arc_lengths[size - 1]) / (arc_lengths[size - 2] - arc_lengths[size - 1])
            };
            const value_type val{lerp(anchor_points[size - 1], anchor_points[size - 2], ratio)};
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return val + normal * l;
        }
        const param_type ratio{
            (s - arc_lengths[pos - 1]) / (arc_lengths[pos] - arc_lengths[pos - 1])
        };
        const auto [val, diff, diff2] = process(pos, ratio);
        const value_type normal{
            diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
        };
        return val + normal * l;
    }

    [[using gnu: pure, always_inline]]
    auto eval(param_type s, param_type l, param_type v) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{arc_lengths.size()};
        const size_type pos = nearestUpperElement(arc_lengths, s) - arc_lengths.begin();
        if (pos < 2) {
            const param_type ratio{(s - arc_lengths[0]) / (arc_lengths[1] - arc_lengths[0])};
            const value_type val{lerp(anchor_points[0], anchor_points[1], ratio)};
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return val + normal * l + binormal * v;
        }
        if (pos == size) {
            const param_type ratio{
                (s - arc_lengths[size - 1]) / (arc_lengths[size - 2] - arc_lengths[size - 1])
            };
            const value_type val{lerp(anchor_points[size - 1], anchor_points[size - 2], ratio)};
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return val + normal * l + binormal * v;
        }
        const param_type ratio{
            (s - arc_lengths[pos - 1]) / (arc_lengths[pos] - arc_lengths[pos - 1])
        };
        const auto [val, diff, diff2] = process(pos, ratio);
        const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
        const value_type binormal{diff.cross(diff2).normalized()};
        return val + normal * l + binormal * v;
    }

    [[using gnu: pure, always_inline]]
    auto eval(SlDuplet<param_type> sl) const noexcept -> value_type
        requires Vec2Arithmetic<value_type>
    {
        return eval(sl.s, sl.l);
    }

    [[using gnu: pure, always_inline]]
    auto eval(SlvTriplet<param_type> slv) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        return eval(slv.s, slv.l, slv.v);
    }

    [[using gnu: pure, always_inline]]
    auto tangent(param_type s) const noexcept -> value_type {
        return m_vec_of_s.derivative(s).normalized();
    }

    [[using gnu: pure, always_inline]]
    auto normal(param_type s) const noexcept -> value_type
        requires Vec2Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{arc_lengths.size()};
        const size_type pos = nearestUpperElement(arc_lengths, s) - arc_lengths.begin();
        if (pos < 2) {
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            return diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1);
        }
        if (pos == size) {
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            return diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1);
        }
        const param_type ratio{
            (s - arc_lengths[pos - 1]) / (arc_lengths[pos] - arc_lengths[pos - 1])
        };
        const auto [val, diff, diff2] = process(pos, ratio);
        return diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1);
    }

    [[using gnu: pure, always_inline]]
    auto normal(param_type s) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{arc_lengths.size()};
        const size_type pos = nearestUpperElement(arc_lengths, s) - arc_lengths.begin();
        if (pos < 2) {
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            return diff.cross(diff2.cross(diff)).normalized();
        }
        if (pos == size) {
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            return diff.cross(diff2.cross(diff)).normalized();
        }
        const auto [val, diff, diff2] = process(pos, 0.0);
        return diff.cross(diff2.cross(diff)).normalized();
    }

    [[using gnu: pure, always_inline]]
    auto binormal(param_type s) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{arc_lengths.size()};
        const size_type pos = nearestUpperElement(arc_lengths, s) - arc_lengths.begin();
        if (pos < 2) {
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            return diff.cross(diff2).normalized();
        }
        if (pos == size) {
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            return diff.cross(diff2).normalized();
        }
        const auto [val, diff, diff2] = process(pos, 0.0);
        return diff.cross(diff2).normalized();
    }

    [[using gnu: pure, always_inline]]
    auto curvature([[maybe_unused]] param_type s) const noexcept -> param_type {
        return std::numeric_limits<param_type>::quiet_NaN();
    }

    [[using gnu: pure, always_inline]]
    auto torsion([[maybe_unused]] param_type s) const noexcept -> param_type
        requires Vec3Arithmetic<value_type>
    {
        return std::numeric_limits<param_type>::quiet_NaN();
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(value_type point) const noexcept -> SlDuplet<param_type>
        requires Vec2Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{anchorPoints().size()};
        const size_type pos = nearestUpperElement(anchor_points, point) - anchor_points.begin();
        if (pos < 2) {
            const value_type r{point - anchor_points[0]};
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return SlDuplet<param_type>{.s{arc_lengths[0] + r.dot(tangent)}, .l{r.dot(normal)}};
        };
        if (pos == size) {
            const value_type r{point - anchor_points[size - 1]};
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return SlDuplet<param_type>{
                .s{arc_lengths[size - 1] + r.dot(tangent)}, .l{r.dot(normal)}
            };
        }
        const value_type r{point - anchor_points[pos - 1]};
        const auto [val, diff, diff2] = process(pos, 0.0);
        const value_type tangent{diff.normalized()};
        const value_type normal{
            diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
        };
        return SlDuplet<param_type>{.s{arc_lengths[pos - 1] + r.dot(tangent)}, .l{r.dot(normal)}};
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(value_type point, param_type start_s, param_type end_s) const noexcept
        -> SlDuplet<param_type>
        requires Vec2Arithmetic<value_type>
    {
        if (start_s > end_s) {
            std::swap(start_s, end_s);
        }
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{anchorPoints().size()};
        const size_type istart =
            nearestUpperElement(arc_lengths, start_s) - arc_lengths.begin() - 1;
        const size_type iend = nearestUpperElement(arc_lengths, end_s) - arc_lengths.begin();
        const size_type pos =
            nearestUpperElement(
                anchor_points.begin() + istart, anchor_points.begin() + iend, point
            ) -
            anchor_points.begin();
        if (pos < istart + 2 && size > istart + 2) {
            const value_type r{point - anchor_points[istart]};
            const value_type diff{anchor_points[istart + 1] - anchor_points[istart]};
            const value_type diff2{
                (anchor_points[istart + 2] - anchor_points[istart + 1]) /
                    (arc_lengths[istart + 2] - arc_lengths[istart + 1]) -
                (anchor_points[istart + 1] - anchor_points[istart]) /
                    (arc_lengths[istart + 1] - arc_lengths[istart])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return SlDuplet<param_type>{
                .s{arc_lengths[istart] + r.dot(tangent)}, .l{r.dot(normal)}
            };
        };
        if (pos == iend) {
            const value_type r{point - anchor_points[iend - 1]};
            const value_type diff{anchor_points[iend - 1] - anchor_points[iend - 2]};
            const value_type diff2{
                (anchor_points[iend - 1] - anchor_points[iend - 2]) /
                    (arc_lengths[iend - 1] - arc_lengths[iend - 2]) -
                (anchor_points[iend - 2] - anchor_points[iend - 3]) /
                    (arc_lengths[iend - 2] - arc_lengths[iend - 3])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{
                diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
            };
            return SlDuplet<param_type>{
                .s{arc_lengths[iend - 1] + r.dot(tangent)}, .l{r.dot(normal)}
            };
        }
        const value_type r{point - anchor_points[pos - 1]};
        const auto [val, diff, diff2] = process(pos, 0.0);
        const value_type tangent{diff.normalized()};
        const value_type normal{
            diff.rotateHalfPi().normalized() * (diff.crossProj(diff2) > 0 ? 1 : -1)
        };
        return SlDuplet<param_type>{.s{arc_lengths[pos - 1] + r.dot(tangent)}, .l{r.dot(normal)}};
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse([[maybe_unused]] value_type point) const noexcept -> SlvTriplet<param_type>
        requires Vec3Arithmetic<value_type>
    {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{anchorPoints().size()};
        const size_type pos = nearestUpperElement(anchor_points, point) - anchor_points.begin();
        if (pos < 2) {
            const value_type r{point - anchor_points[0]};
            const value_type diff{anchor_points[1] - anchor_points[0]};
            const value_type diff2{
                (anchor_points[2] - anchor_points[1]) / (arc_lengths[2] - arc_lengths[1]) -
                (anchor_points[1] - anchor_points[0]) / (arc_lengths[1] - arc_lengths[0])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return SlvTriplet<param_type>{
                .s{arc_lengths[0] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
            };
        }
        if (pos == size) {
            const value_type r{point - anchor_points[size - 1]};
            const value_type diff{anchor_points[size - 1] - anchor_points[size - 2]};
            const value_type diff2{
                (anchor_points[size - 1] - anchor_points[size - 2]) /
                    (arc_lengths[size - 1] - arc_lengths[size - 2]) -
                (anchor_points[size - 2] - anchor_points[size - 3]) /
                    (arc_lengths[size - 2] - arc_lengths[size - 3])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return SlvTriplet<param_type>{
                .s{arc_lengths[size - 1] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
            };
        }
        const value_type r{point - anchor_points[pos - 1]};
        const auto [val, diff, diff2] = process(pos, 0.0);
        const value_type tangent{diff.normalized()};
        const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
        const value_type binormal{diff.cross(diff2).normalized()};
        return SlvTriplet<param_type>{
            .s{arc_lengths[pos - 1] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
        };
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(value_type point, param_type start_s, param_type end_s) const noexcept
        -> SlvTriplet<param_type>
        requires Vec3Arithmetic<value_type>
    {
        if (start_s > end_s) {
            std::swap(start_s, end_s);
        }
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const size_type size{anchorPoints().size()};
        const size_type istart =
            nearestUpperElement(arc_lengths, start_s) - arc_lengths.begin() - 1;
        const size_type iend = nearestUpperElement(arc_lengths, end_s) - arc_lengths.begin();
        const size_type pos =
            nearestUpperElement(
                anchor_points.begin() + istart, anchor_points.begin() + iend, point
            ) -
            anchor_points.begin();
        if (pos < istart + 2 && size > istart + 2) {
            const value_type r{point - anchor_points[istart]};
            const value_type diff{anchor_points[istart + 1] - anchor_points[istart]};
            const value_type diff2{
                (anchor_points[istart + 2] - anchor_points[istart + 1]) /
                    (arc_lengths[istart + 2] - arc_lengths[istart + 1]) -
                (anchor_points[istart + 1] - anchor_points[istart]) /
                    (arc_lengths[istart + 1] - arc_lengths[istart])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return SlvTriplet<param_type>{
                .s{arc_lengths[istart] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
            };
        }
        if (pos == iend) {
            const value_type r{point - anchor_points[iend - 1]};
            const value_type diff{anchor_points[iend - 1] - anchor_points[iend - 2]};
            const value_type diff2{
                (anchor_points[iend - 1] - anchor_points[iend - 2]) /
                    (arc_lengths[iend - 1] - arc_lengths[iend - 2]) -
                (anchor_points[iend - 2] - anchor_points[iend - 3]) /
                    (arc_lengths[iend - 2] - arc_lengths[iend - 3])
            };
            const value_type tangent{diff.normalized()};
            const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
            const value_type binormal{diff.cross(diff2).normalized()};
            return SlvTriplet<param_type>{
                .s{arc_lengths[iend - 1] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
            };
        }
        const value_type r{point - anchor_points[pos - 1]};
        const auto [val, diff, diff2] = process(pos, 0.0);
        const value_type tangent{diff.normalized()};
        const value_type normal{diff.cross(diff2.cross(diff)).normalized()};
        const value_type binormal{diff.cross(diff2).normalized()};
        return SlvTriplet<param_type>{
            .s{arc_lengths[pos - 1] + r.dot(tangent)}, .l{r.dot(normal)}, .v{r.dot(binormal)}
        };
    }

    [[using gnu: pure, always_inline]]
    auto operator()(param_type s) const noexcept -> value_type {
        return eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(param_type s, param_type l) const noexcept -> value_type
        requires Vec2Arithmetic<value_type>
    {
        return eval(s, l);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(param_type s, param_type l, param_type v) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        return eval(s, l, v);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(SlDuplet<param_type> sl) const noexcept -> value_type
        requires Vec2Arithmetic<value_type>
    {
        return eval(sl);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(SlvTriplet<param_type> slv) const noexcept -> value_type
        requires Vec3Arithmetic<value_type>
    {
        return eval(slv);
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
    auto arcLengths() const noexcept -> std::span<const param_type> {
        return m_vec_of_s.ts();
    }

    [[using gnu: pure, always_inline]]
    auto anchorPoints() const noexcept -> std::span<const value_type> {
        return m_vec_of_s.ys();
    }

  private:
    using param_allocator_type =
        typename std::allocator_traits<allocator_type>::template rebind_alloc<param_type>;

    [[using gnu: pure, flatten, leaf, hot]]
    auto process(size_type pos, param_type ratio) const noexcept
        -> std::tuple<value_type, value_type, value_type> {
        auto arc_lengths{arcLengths()};
        auto anchor_points{anchorPoints()};
        const value_type val{lerp(anchor_points[pos - 1], anchor_points[pos], ratio)};
        const value_type diff{anchor_points[pos] - anchor_points[pos - 1]};
        const value_type diff2{
            (anchor_points[pos] - anchor_points[pos - 1]) /
                (arc_lengths[pos] - arc_lengths[pos - 1]) -
            (anchor_points[pos - 1] - anchor_points[pos - 2]) /
                (arc_lengths[pos - 1] - arc_lengths[pos - 2])
        };
        return {val, diff, diff2};
    }

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_vec_of_s;
        return;
    }

    PiecewiseLinearFunction<value_type, allocator_type> m_vec_of_s;
};

template <Allocatory Alloc = ::boyle::common::AlignedAllocator<Vec2s, alignof(Vec2s)>>
using PiecewiseLinearCurve2s = PiecewiseLinearCurve<Vec2s, Alloc>;

template <Allocatory Alloc = ::boyle::common::AlignedAllocator<Vec2d, alignof(Vec2d)>>
using PiecewiseLinearCurve2d = PiecewiseLinearCurve<Vec2d, Alloc>;

template <Allocatory Alloc = ::boyle::common::AlignedAllocator<Vec3s, alignof(Vec3s)>>
using PiecewiseLinearCurve3s = PiecewiseLinearCurve<Vec3s, Alloc>;

template <Allocatory Alloc = ::boyle::common::AlignedAllocator<Vec3d, alignof(Vec3d)>>
using PiecewiseLinearCurve3d = PiecewiseLinearCurve<Vec3d, Alloc>;

namespace pmr {

template <VecArithmetic T>
using PiecewiseLinearCurve =
    boyle::math::PiecewiseLinearCurve<T, std::pmr::polymorphic_allocator<T>>;

using PiecewiseLinearCurve2s =
    ::boyle::math::PiecewiseLinearCurve<Vec2s, std::pmr::polymorphic_allocator<Vec2s>>;

using PiecewiseLinearCurve2d =
    ::boyle::math::PiecewiseLinearCurve<Vec2d, std::pmr::polymorphic_allocator<Vec2d>>;

using PiecewiseLinearCurve3s =
    ::boyle::math::PiecewiseLinearCurve<Vec3s, std::pmr::polymorphic_allocator<Vec3s>>;

using PiecewiseLinearCurve3d =
    ::boyle::math::PiecewiseLinearCurve<Vec3d, std::pmr::polymorphic_allocator<Vec3d>>;

} // namespace pmr

} // namespace boyle::math
