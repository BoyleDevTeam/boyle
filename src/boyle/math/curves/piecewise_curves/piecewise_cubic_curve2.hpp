/**
 * @file piecewise_cubic_curve2.hpp
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
#include <array>
#include <concepts>
#include <format>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/cubic_interpolation.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_cubic_function1.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::math {

template <std::floating_point T>
class [[nodiscard]] PiecewiseCubicCurve2 final {
    friend class boost::serialization::access;

  public:
    using value_type = Vec2<T>;
    using param_type = T;
    using BoundaryMode = typename PiecewiseCubicFunction1<Vec2<T>, T>::BoundaryMode;

    static constexpr T kDuplicateCriterion{kEpsilon};

    [[using gnu: always_inline]]
    explicit PiecewiseCubicCurve2(std::vector<Vec2<T>> anchor_points, T s0 = 0.0)
        : PiecewiseCubicCurve2(
              std::move(anchor_points), BoundaryMode{2, Vec2<T>{0.0, 0.0}},
              BoundaryMode{2, Vec2<T>{0.0, 0.0}}, s0
          ) {}

    [[using gnu: flatten, leaf]]
    explicit PiecewiseCubicCurve2(
        std::vector<Vec2<T>> anchor_points, BoundaryMode b0, BoundaryMode bf, T s0 = 0.0
    ) noexcept(!BOYLE_CHECK_PARAMS) {
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
        std::vector<T> temp_arc_lengths;
        temp_arc_lengths.reserve(size);
        temp_arc_lengths.push_back(s0);
        for (typename std::vector<Vec2<T>>::const_iterator it = anchor_points.cbegin() + 1;
             it != anchor_points.cend(); ++it) {
            temp_arc_lengths.push_back(temp_arc_lengths.back() + (it - 1)->euclideanTo(*it));
        }
        const PiecewiseCubicFunction1<Vec2<T>, T> temp_vec2_of_s{
            temp_arc_lengths, anchor_points, b0, bf
        };

        const std::vector<Vec2<T>>& ddys = temp_vec2_of_s.ddys();
        std::vector<T> arc_lengths;
        arc_lengths.reserve(size);
        arc_lengths.push_back(s0);
        for (std::size_t i = 1; i < size; ++i) {
            arc_lengths.push_back(
                arc_lengths.back() + calcArcLength(
                                         anchor_points[i - 1], anchor_points[i], ddys[i - 1],
                                         ddys[i], temp_arc_lengths[i] - temp_arc_lengths[i - 1]
                                     )
            );
        }
        m_vec2_of_s = PiecewiseCubicFunction1<Vec2<T>, T>{
            std::move(arc_lengths), std::move(anchor_points), b0, bf
        };
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseCubicCurve2);
    ~PiecewiseCubicCurve2() noexcept = default;

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

    [[using gnu: pure]]
    auto curvature(T s) const noexcept -> T {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
        const std::size_t size = arc_lengths.size();
        const std::size_t pos =
            nearestUpperElement(arc_lengths.cbegin(), arc_lengths.cend(), s) - arc_lengths.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! s should be greater than minS(): s = {0:.6f} while "
                "minS() = {1:.6f}. Use extra interpolating!",
                s, arc_lengths.front()
            );
            return 0.0;
        }
        if (pos == size) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! s should be less than maxS(): s = {0:.6f} while "
                "maxS() = {1:.6f}. Use extra interpolating!",
                s, arc_lengths.back()
            );
            return 0.0;
        }
        const T h = arc_lengths[pos] - arc_lengths[pos - 1];
        const T ratio = (s - arc_lengths[pos - 1]) / h;
        const Vec2<T> derivative =
            cuberpd(anchor_points[pos - 1], anchor_points[pos], ddys[pos - 1], ddys[pos], ratio, h);
        const Vec2<T> derivative2 = lerp(ddys[pos - 1], ddys[pos], ratio);
        const T derivative_norm = derivative.norm();
        return derivative.cross(derivative2) /
               (derivative_norm * derivative_norm * derivative_norm);
    }

    [[using gnu: pure, flatten, leaf]]
    auto inverse(Vec2<T> point) const noexcept -> SlDuplet<T> {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
        const std::size_t pos =
            nearestUpperElement(anchor_points.cbegin(), anchor_points.cend(), point) -
            anchor_points.cbegin();
        Vec2<T> derivative, derivative2, r;
        T h, derivative_norm, s, l;
        if (pos == 0) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = cuberpd(anchor_points[0], anchor_points[1], ddys[0], ddys[1], 0.0, h);
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else if (pos == anchor_points.size()) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = cuberpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else {
            const Vec2<T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<T>& upper_derivative2{ddys[pos]};
            T ratio = 0.5;
            h = m_vec2_of_s.ts()[pos] - m_vec2_of_s.ts()[pos - 1];
            r = point - cuberp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, ratio, h
                        );
            derivative = cuberpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2, ratio,
                h
            );
            derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            for (std::size_t num_iter = 3; num_iter != 0U; --num_iter) {
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
            l = derivative.cross(r) / derivative.norm();
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
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
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
        Vec2<T> derivative, derivative2, r;
        T h, derivative_norm, s, l;
        if (pos == istart) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = cuberpd(anchor_points[0], anchor_points[1], ddys[0], ddys[1], 0.0, h);
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else if (pos == iend) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = cuberpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else {
            const Vec2<T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<T>& upper_derivative2{ddys[pos]};
            T ratio = 0.5;
            h = m_vec2_of_s.ts()[pos] - m_vec2_of_s.ts()[pos - 1];
            r = point - cuberp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, ratio, h
                        );
            derivative = cuberpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2, ratio,
                h
            );
            derivative2 = lerp(lower_derivative2, upper_derivative2, ratio);
            for (std::size_t num_iter = 3; num_iter != 0U; --num_iter) {
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
            l = derivative.cross(r) / derivative.norm();
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
    [[using gnu: const, always_inline]]
    static auto calcArcLength(
        Vec2<T> start, Vec2<T> end, Vec2<T> ddstart, Vec2<T> ddend, T scale
    ) noexcept -> T {
        /* using boost::math::ccmath::sqrt;
        constexpr std::size_t kNumGaussLegendreKnots{5};
        constexpr std::array<T, 5> kGaussLegendreKnots{
            0.5 - sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 6.0,
            0.5 - sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 6.0, 0.5,
            0.5 + sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 6.0,
            0.5 + sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 6.0};
        constexpr std::array<T, 5> kGaussLegendreWeights{
            (322.0 - 13.0 * sqrt(70.0)) / 1800.0, (322.0 + 13.0 * sqrt(70.0)) / 1800.0,
            128.0 / 450.0, (322.0 + 13.0 * sqrt(70.0)) / 1800.0,
            (322.0 - 13.0 * sqrt(70.0)) / 1800.0}; */
        constexpr std::size_t kNumGaussLegendreKnots{5};
        constexpr std::array<T, 5> kGaussLegendreKnots{
            0.04691007703066800360, 0.23076534494715845448, 0.5, 0.76923465505284154551,
            0.95308992296933199639
        };
        constexpr std::array<T, 5> kGaussLegendreWeights{
            0.11846344252809454375, 0.23931433524968323402, 128.0 / 450.0, 0.23931433524968323402,
            0.11846344252809454375
        };
        T result{0.0};
        for (std::size_t i = 0; i < kNumGaussLegendreKnots; ++i) {
            result += cuberpd(start, end, ddstart, ddend, kGaussLegendreKnots[i], scale).norm() *
                      kGaussLegendreWeights[i];
        }
        result *= scale;
        return result;
    }

    [[using gnu: always_inline]]
    explicit PiecewiseCubicCurve2(PiecewiseCubicFunction1<math::Vec2<T>, T> vec2_of_s) noexcept
        : m_vec2_of_s{std::move(vec2_of_s)} {}

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_vec2_of_s;
        return;
    }

    PiecewiseCubicFunction1<Vec2<T>, T> m_vec2_of_s{};
};

using PiecewiseCubicCurve2f = PiecewiseCubicCurve2<float>;
using PiecewiseCubicCurve2d = PiecewiseCubicCurve2<double>;

} // namespace boyle::math
