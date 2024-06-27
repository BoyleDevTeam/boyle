/**
 * @file piecewise_linear_curve2.hpp
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
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_linear_function1.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::math {

template <std::floating_point T>
class [[nodiscard]] PiecewiseLinearCurve2 final {
    friend class boost::serialization::access;

  public:
    using value_type = Vec2<T>;
    using param_type = T;

    static constexpr T kDuplicateCriterion{kEpsilon};

    [[using gnu: flatten, leaf]]
    explicit PiecewiseLinearCurve2(std::vector<Vec2<T>> anchor_points, T s0 = 0.0) noexcept(
        !BOYLE_CHECK_PARAMS
    ) {
#if BOYLE_CHECK_PARAMS == 1
        if (anchor_points.size() < 2) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! sizes of anchor_points must be greater than 2: "
                "anchor_points.size() = {0:d}",
                anchor_points.size()
            ));
        }
#endif
        const std::size_t size = anchor_points.size();
        std::vector<T> arc_lengths;
        arc_lengths.reserve(size);
        arc_lengths.push_back(s0);
        for (typename std::vector<Vec2<T>>::const_iterator it = anchor_points.cbegin() + 1;
             it != anchor_points.cend(); ++it) {
            arc_lengths.push_back(arc_lengths.back() + (it - 1)->euclideanTo(*it));
        }
        m_vec2_of_s =
            PiecewiseLinearFunction1<Vec2<T>, T>{std::move(arc_lengths), std::move(anchor_points)};
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseLinearCurve2);
    ~PiecewiseLinearCurve2() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto eval(T s) const noexcept -> Vec2<T> {
        return m_vec2_of_s.eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto eval(T s, T l) const noexcept -> Vec2<T> {
        return eval(s) + l * orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    auto eval(SlDuplet<T> sl) const noexcept -> Vec2<T> {
        return eval(sl.s) + sl.l * orthonormal(sl.s);
    }

    [[using gnu: pure, always_inline]]
    auto tangent(T s) const noexcept -> Vec2<T> {
        return normalize(m_vec2_of_s.derivative(s));
    }

    [[using gnu: pure, always_inline]]
    auto orthonormal(T s) const noexcept -> Vec2<T> {
        return tangent(s).selfRotateHalfPi();
    }

    [[using gnu: pure, flatten, leaf]]
    auto inverse(Vec2<T> point) const noexcept -> SlDuplet<T> {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::size_t pos =
            nearestUpperElement(anchor_points.cbegin(), anchor_points.cend(), point) -
            anchor_points.cbegin();
        Vec2<T> diff, r;
        T diff_norm, s, l;
        if (pos == 0) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_norm = diff.norm();
            s = arc_lengths.front() + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        } else if (pos == anchor_points.size()) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_norm = diff.norm();
            s = arc_lengths.back() + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_norm = diff.norm();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        }
        return SlDuplet<T>{s, l};
    }

    [[using gnu: pure]]
    auto inverse(Vec2<T> point, T start_s, T end_s) const noexcept -> SlDuplet<T> {
        if (start_s > end_s) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The start_s should always be less than end_s: "
                "start_s = {0:.6f} while end_s = {1:.6f}.",
                start_s, end_s
            );
            std::swap(start_s, end_s);
        }
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::size_t istart =
            nearestUpperElement(arc_lengths.cbegin(), arc_lengths.cend(), start_s) -
            arc_lengths.cbegin() - 1;
        const std::size_t iend =
            nearestUpperElement(arc_lengths.cbegin(), arc_lengths.cend(), end_s) -
            arc_lengths.cbegin();
        const std::size_t pos =
            nearestUpperElement(
                anchor_points.cbegin() + istart, anchor_points.cend() + iend, point
            ) -
            anchor_points.cbegin();
        Vec2<T> diff, r;
        T diff_norm, s, l;
        if (pos == istart) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_norm = diff.norm();
            s = arc_lengths.front() + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        } else if (pos == iend) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_norm = diff.norm();
            s = arc_lengths.back() + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_norm = diff.norm();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_norm;
            l = diff.cross(r) / diff_norm;
        }
        return SlDuplet<T>{s, l};
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto eval(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<Vec2<T>> points;
        points.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            points.emplace_back(eval(*it));
        }
        return points;
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto eval(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<Vec2<T>>
        requires std::same_as<SlDuplet<T>, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<Vec2<T>> points;
        points.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            points.emplace_back(eval(*it));
        }
        return points;
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto tangent(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<Vec2<T>> tangents;
        tangents.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            tangents.emplace_back(tangent(*it));
        }
        return tangents;
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto orthonormal(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<Vec2<T>>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<Vec2<T>> orthonormals;
        orthonormals.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            orthonormals.emplace_back(orthonormal(*it));
        }
        return orthonormals;
    }

    [[using gnu: pure, always_inline]]
    auto operator()(T s) const noexcept -> Vec2<T> {
        return eval(s);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(T s, T l) const noexcept -> Vec2<T> {
        return eval(s, l);
    }

    [[using gnu: pure, always_inline]]
    auto operator()(SlDuplet<T> sl) const noexcept -> Vec2<T> {
        return eval(sl);
    }

    [[using gnu: pure, always_inline]]
    auto minS() const noexcept -> T {
        return m_vec2_of_s.minT();
    }

    [[using gnu: pure, always_inline]]
    auto maxS() const noexcept -> T {
        return m_vec2_of_s.maxT();
    }

    [[using gnu: pure, always_inline]]
    auto front() const noexcept -> Vec2<T> {
        return m_vec2_of_s.ys().front();
    }

    [[using gnu: pure, always_inline]]
    auto back() const noexcept -> Vec2<T> {
        return m_vec2_of_s.ys().back();
    }

    [[using gnu: pure, always_inline]]
    auto arcLengths() const noexcept -> const std::vector<T>& {
        return m_vec2_of_s.ts();
    }

    [[using gnu: pure, always_inline]]
    auto anchorPoints() const noexcept -> const std::vector<Vec2<T>>& {
        return m_vec2_of_s.ys();
    }

  private:
    [[using gnu: always_inline]]
    explicit PiecewiseLinearCurve2(PiecewiseLinearFunction1<math::Vec2<T>, T> vec2_of_s)
        : m_vec2_of_s{std::move(vec2_of_s)} {}

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_vec2_of_s;
        return;
    }

    PiecewiseLinearFunction1<Vec2<T>, T> m_vec2_of_s{};
};

using PiecewiseLinearCurve2f = PiecewiseLinearCurve2<float>;
using PiecewiseLinearCurve2d = PiecewiseLinearCurve2<double>;

} // namespace boyle::math
