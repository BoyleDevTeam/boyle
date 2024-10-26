/**
 * @file piecewise_quintic_curve2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-03
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "fmt/format.h"

#include "boyle/math/cubic_interpolation.hpp"
#include "boyle/math/duplet.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_quintic_function1.hpp"
#include "boyle/math/quintic_interpolation.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::math {

template <std::floating_point T>
class [[nodiscard]] PiecewiseQuinticCurve2 final {
    friend class boost::serialization::access;

  public:
    using value_type = Vec2<T>;
    using param_type = T;
    using BoundaryMode = typename PiecewiseQuinticFunction1<Vec2<T>, T>::BoundaryMode;

    static constexpr T kDuplicateCriterion{1E-8};

    PiecewiseQuinticCurve2() noexcept = default;
    PiecewiseQuinticCurve2(const PiecewiseQuinticCurve2& other) noexcept = default;
    auto operator=(const PiecewiseQuinticCurve2& other
    ) noexcept -> PiecewiseQuinticCurve2& = default;
    PiecewiseQuinticCurve2(PiecewiseQuinticCurve2&& other) noexcept = default;
    auto operator=(PiecewiseQuinticCurve2&& other) noexcept -> PiecewiseQuinticCurve2& = default;
    ~PiecewiseQuinticCurve2() noexcept = default;

    [[using gnu: always_inline]]
    explicit PiecewiseQuinticCurve2(std::vector<Vec2<T>> anchor_points, T s0 = 0.0)
        : PiecewiseQuinticCurve2(
              std::move(anchor_points),
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, Vec2<T>{0.0, 0.0}}, BoundaryMode{4, Vec2<T>{0.0, 0.0}}
              },
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, Vec2<T>{0.0, 0.0}}, BoundaryMode{4, Vec2<T>{0.0, 0.0}}
              },
              s0
          ) {}

    [[using gnu: ]]
    explicit PiecewiseQuinticCurve2(
        std::vector<Vec2<T>> anchor_points, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf, T s0 = 0.0
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
        std::vector<T> temp_arc_lengths;
        temp_arc_lengths.reserve(size);
        temp_arc_lengths.push_back(s0);
        for (typename std::vector<Vec2<T>>::const_iterator it = anchor_points.cbegin() + 1;
             it != anchor_points.cend(); ++it) {
            temp_arc_lengths.push_back(temp_arc_lengths.back() + (it - 1)->euclideanTo(*it));
        }
        const PiecewiseQuinticFunction1<Vec2<T>, T> temp_vec2_of_s{
            temp_arc_lengths, anchor_points, b0, bf
        };

        const std::vector<Vec2<T>>& ddys = temp_vec2_of_s.ddys();
        const std::vector<Vec2<T>>& d4ys = temp_vec2_of_s.d4ys();
        std::vector<T> arc_lengths;
        arc_lengths.reserve(size);
        arc_lengths.push_back(s0);
        for (std::size_t i = 1; i < size; ++i) {
            arc_lengths.push_back(
                arc_lengths.back() + ::boyle::math::calcArcLength(
                                         anchor_points[i - 1], anchor_points[i], ddys[i - 1],
                                         ddys[i], d4ys[i - 1], d4ys[i],
                                         temp_arc_lengths[i] - temp_arc_lengths[i - 1]
                                     )
            );
        }
        m_vec2_of_s = PiecewiseQuinticFunction1<Vec2<T>, T>{
            std::move(arc_lengths), std::move(anchor_points), b0, bf
        };
    }

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

    [[using gnu: pure, flatten, leaf, hot]]
    auto curvature(T s) const noexcept -> T {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
        const std::vector<Vec2<T>>& d4ys = m_vec2_of_s.d4ys();
        const std::size_t size = arc_lengths.size();
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
        const T h = arc_lengths[pos] - arc_lengths[pos - 1];
        const T ratio = (s - arc_lengths[pos - 1]) / h;
        const Vec2<T> derivative = quinerpd(
            anchor_points[pos - 1], anchor_points[pos], ddys[pos - 1], ddys[pos], d4ys[pos - 1],
            d4ys[pos - 1], ratio, h
        );
        const Vec2<T> derivative2 =
            cuberp(ddys[pos - 1], ddys[pos], d4ys[pos - 1], d4ys[pos], ratio, h);
        const T derivative_norm = derivative.norm();
        return derivative.cross(derivative2) /
               (derivative_norm * derivative_norm * derivative_norm);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(Vec2<T> point) const noexcept -> SlDuplet<T> {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
        const std::vector<Vec2<T>>& d4ys = m_vec2_of_s.d4ys();
        const std::size_t pos =
            nearestUpperElement(
                std::ranges::subrange{anchor_points.cbegin(), anchor_points.cend()}, point
            ) -
            anchor_points.cbegin();
        Vec2<T> derivative, derivative2, r;
        T h, derivative_norm, s, l;
        if (pos == 0) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = quinerpd(
                anchor_points[0], anchor_points[1], ddys[0], ddys[1], d4ys[0], d4ys[1], 0.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else if (pos == anchor_points.size()) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = quinerpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1],
                d4ys[pos - 2], d4ys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else {
            const Vec2<T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<T>& upper_derivative2{ddys[pos]};
            const Vec2<T>& lower_derivative4{d4ys[pos - 1]};
            const Vec2<T>& upper_derivative4{d4ys[pos]};
            T ratio = 0.5;
            h = m_vec2_of_s.ts()[pos] - m_vec2_of_s.ts()[pos - 1];
            r = point - quinerp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
                        );
            derivative = quinerpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                lower_derivative4, upper_derivative4, ratio, h
            );
            derivative2 = cuberp(
                lower_derivative2, upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
            );
            for (std::size_t num_iter = 3; num_iter != 0U; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.normSqr()) * h);
                r = point - quinerp(
                                lower_anchor_point, upper_anchor_point, lower_derivative2,
                                upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
                            );
                derivative = quinerpd(
                    lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                    lower_derivative4, upper_derivative4, ratio, h
                );
                derivative2 = cuberp(
                    lower_derivative2, upper_derivative2, lower_derivative4, upper_derivative4,
                    ratio, h
                );
            }
            s = lerp(arc_lengths[pos - 1], arc_lengths[pos], ratio);
            l = derivative.cross(r) / derivative.norm();
        }
        return SlDuplet<T>{s, l};
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto inverse(Vec2<T> point, T start_s, T end_s) const noexcept -> SlDuplet<T> {
        if (start_s > end_s) {
            std::swap(start_s, end_s);
        }
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<T>>& ddys = m_vec2_of_s.ddys();
        const std::vector<Vec2<T>>& d4ys = m_vec2_of_s.d4ys();
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
                std::ranges::subrange{anchor_points.cbegin() + istart, anchor_points.cend() + iend},
                point
            ) -
            anchor_points.cbegin();
        Vec2<T> derivative, derivative2, r;
        T h, derivative_norm, s, l;
        if (pos == istart) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = quinerpd(
                anchor_points[0], anchor_points[1], ddys[0], ddys[1], d4ys[0], d4ys[1], 0.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.front() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else if (pos == iend) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = quinerpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1],
                d4ys[pos - 2], d4ys[pos - 1], 1.0, h
            );
            derivative_norm = derivative.norm();
            s = arc_lengths.back() + derivative.dot(r) / derivative_norm;
            l = derivative.cross(r) / derivative_norm;
        } else {
            const Vec2<T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<T>& upper_derivative2{ddys[pos]};
            const Vec2<T>& lower_derivative4{d4ys[pos - 1]};
            const Vec2<T>& upper_derivative4{d4ys[pos]};
            T ratio = 0.5;
            h = m_vec2_of_s.ts()[pos] - m_vec2_of_s.ts()[pos - 1];
            r = point - quinerp(
                            lower_anchor_point, upper_anchor_point, lower_derivative2,
                            upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
                        );
            derivative = quinerpd(
                lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                lower_derivative4, upper_derivative4, ratio, h
            );
            derivative2 = cuberp(
                lower_derivative2, upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
            );
            for (std::size_t num_iter = 3; num_iter != 0U; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.normSqr()) * h);
                r = point - quinerp(
                                lower_anchor_point, upper_anchor_point, lower_derivative2,
                                upper_derivative2, lower_derivative4, upper_derivative4, ratio, h
                            );
                derivative = quinerpd(
                    lower_anchor_point, upper_anchor_point, lower_derivative2, upper_derivative2,
                    lower_derivative4, upper_derivative4, ratio, h
                );
                derivative2 = cuberp(
                    lower_derivative2, upper_derivative2, lower_derivative4, upper_derivative4,
                    ratio, h
                );
            }
            s = lerp(arc_lengths[pos - 1], arc_lengths[pos], ratio);
            l = derivative.cross(r) / derivative.norm();
        }
        return SlDuplet<T>{s, l};
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
    explicit PiecewiseQuinticCurve2(PiecewiseQuinticFunction1<math::Vec2<T>, T> vec2_of_s) noexcept
        : m_vec2_of_s{std::move(vec2_of_s)} {}

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_vec2_of_s;
        return;
    }

    PiecewiseQuinticFunction1<Vec2<T>, T> m_vec2_of_s{};
};

using PiecewiseQuinticCurve2f = PiecewiseQuinticCurve2<float>;
using PiecewiseQuinticCurve2d = PiecewiseQuinticCurve2<double>;

} // namespace boyle::math
