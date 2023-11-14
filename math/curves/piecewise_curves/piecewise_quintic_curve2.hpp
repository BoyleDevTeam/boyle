/**
 * @file piecewise_quintic_curve2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-03
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <format>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "spdlog/spdlog.h"

#include "common/utils/macros.hpp"
#include "math/cubic_interpolation.hpp"
#include "math/curves/curve2.hpp"
#include "math/functions/piecewise_functions/piecewise_quintic_function1.hpp"
#include "math/quintic_interpolation.hpp"
#include "math/utils.hpp"
#include "math/vec2.hpp"

namespace tiny_pnc {
namespace math {

template <typename T>
class [[nodiscard]] PiecewiseQuinticCurve2 final : public Curve2<PiecewiseQuinticCurve2<T>, T> {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    using BoundaryMode = typename PiecewiseQuinticFunction1<Vec2<Vec2Mode::XY, T>, T>::BoundaryMode;

    static constexpr T kDuplicateCriterion{kEpsilon};

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseQuinticCurve2);

    [[using gnu: always_inline]] explicit PiecewiseQuinticCurve2(
        std::vector<Vec2<Vec2Mode::XY, T>> anchor_points, T s0 = 0.0
    )
        : PiecewiseQuinticCurve2(
              std::move(anchor_points),
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, Vec2<Vec2Mode::XY, T>{0.0, 0.0}},
                  BoundaryMode{4, Vec2<Vec2Mode::XY, T>{0.0, 0.0}}},
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, Vec2<Vec2Mode::XY, T>{0.0, 0.0}},
                  BoundaryMode{4, Vec2<Vec2Mode::XY, T>{0.0, 0.0}}},
              s0
          ) {}

    [[using gnu: flatten, leaf]] explicit PiecewiseQuinticCurve2(
        std::vector<Vec2<Vec2Mode::XY, T>> anchor_points, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf, T s0 = 0.0
    ) {
        if (anchor_points.size() < 2) {
            std::string error_msg = std::format(
                "Invalid arguments detected! sizes of anchor_points must be greater than 2: "
                "anchor_points.size() = {0:d}",
                anchor_points.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        const std::size_t size = anchor_points.size();
        std::vector<T> temp_arc_lengths;
        temp_arc_lengths.reserve(size);
        temp_arc_lengths.push_back(s0);
        for (typename std::vector<Vec2<Vec2Mode::XY, T>>::const_iterator it =
                 anchor_points.cbegin() + 1;
             it != anchor_points.cend(); ++it) {
            temp_arc_lengths.push_back(temp_arc_lengths.back() + (it - 1)->distanceTo(*it));
        }
        PiecewiseQuinticFunction1<Vec2<Vec2Mode::XY, T>, T> temp_vec2_of_s{
            temp_arc_lengths, anchor_points, b0, bf};

        const std::vector<Vec2<Vec2Mode::XY, T>>& ddys = temp_vec2_of_s.ddys();
        const std::vector<Vec2<Vec2Mode::XY, T>>& d4ys = temp_vec2_of_s.d4ys();
        std::vector<T> arc_lengths;
        arc_lengths.reserve(size);
        arc_lengths.push_back(s0);
        for (std::size_t i = 1; i < size; ++i) {
            arc_lengths.push_back(
                arc_lengths.back() + calcArcLength(
                                         anchor_points[i - 1], anchor_points[i], ddys[i - 1],
                                         ddys[i], d4ys[i - 1], d4ys[i],
                                         temp_arc_lengths[i] - temp_arc_lengths[i - 1]
                                     )
            );
        }
        vec2_of_s_ = PiecewiseQuinticFunction1<Vec2<Vec2Mode::XY, T>, T>{
            std::move(arc_lengths), std::move(anchor_points), b0, bf};
    }

    ~PiecewiseQuinticCurve2() noexcept override = default;

    [[using gnu: pure, always_inline]]
    Vec2<Vec2Mode::XY, T>
    operator()(T s) const noexcept {
        return vec2_of_s_(s);
    }

    [[using gnu: pure, always_inline]]
    Vec2<Vec2Mode::XY, T> tangent(T s) const noexcept {
        return normalize(vec2_of_s_.derivative(s));
    }

    [[using gnu: pure]]
    T curvature(T s) const noexcept {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<Vec2Mode::XY, T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<Vec2Mode::XY, T>>& ddys = vec2_of_s_.ddys();
        const std::vector<Vec2<Vec2Mode::XY, T>>& d4ys = vec2_of_s_.d4ys();
        const std::size_t size = arc_lengths.size();
        const std::size_t pos =
            std::upper_bound(arc_lengths.cbegin(), arc_lengths.cend(), s) - arc_lengths.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! s should be greater than minS(): s = {0:.6f} while "
                "minS() = {1:.6f}. Use extra interpolating!",
                s, arc_lengths.front()
            );
            return 0.0;
        } else if (pos == size) {
            spdlog::warn(
                "Out of range issue detected! s should be less than maxS(): s = {0:.6f} while "
                "maxS() = {1:.6f}. Use extra interpolating!",
                s, arc_lengths.back()
            );
            return 0.0;
        } else {
            const T h = arc_lengths[pos] - arc_lengths[pos - 1];
            const T ratio = (s - arc_lengths[pos - 1]) / h;
            const Vec2<Vec2Mode::XY, T> derivative = quinerpd(
                anchor_points[pos - 1], anchor_points[pos], ddys[pos - 1], ddys[pos], d4ys[pos - 1],
                d4ys[pos - 1], ratio, h
            );
            const Vec2<Vec2Mode::XY, T> derivative2 =
                cuberp(ddys[pos - 1], ddys[pos], d4ys[pos - 1], d4ys[pos], ratio, h);
            const T derivative_length = derivative.length();
            return derivative.cross(derivative2) /
                   (derivative_length * derivative_length * derivative_length);
        }
    }

    [[using gnu: pure, flatten, leaf]]
    Vec2<Vec2Mode::SL, T> inverse(Vec2<Vec2Mode::XY, T> point) const noexcept {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<Vec2Mode::XY, T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<Vec2Mode::XY, T>>& ddys = vec2_of_s_.ddys();
        const std::vector<Vec2<Vec2Mode::XY, T>>& d4ys = vec2_of_s_.d4ys();
        const std::size_t pos =
            closetUpperElement(anchor_points.cbegin(), anchor_points.cend(), point) -
            anchor_points.cbegin();
        Vec2<Vec2Mode::XY, T> derivative, derivative2, r;
        T h, derivative_length, s, l;
        if (pos == 0) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = quinerpd(
                anchor_points[0], anchor_points[1], ddys[0], ddys[1], d4ys[0], d4ys[1], 0.0, h
            );
            derivative_length = derivative.length();
            s = arc_lengths.front() + derivative.dot(r) / derivative_length;
            l = derivative.cross(r) / derivative_length;
        } else if (pos == anchor_points.size()) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = quinerpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1],
                d4ys[pos - 2], d4ys[pos - 1], 1.0, h
            );
            derivative_length = derivative.length();
            s = arc_lengths.back() + derivative.dot(r) / derivative_length;
            l = derivative.cross(r) / derivative_length;
        } else {
            const Vec2<Vec2Mode::XY, T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<Vec2Mode::XY, T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_derivative2{ddys[pos]};
            const Vec2<Vec2Mode::XY, T>& lower_derivative4{d4ys[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_derivative4{d4ys[pos]};
            T ratio = 0.5;
            h = vec2_of_s_.ts()[pos] - vec2_of_s_.ts()[pos - 1];
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
            for (std::size_t num_iter = 3; num_iter; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.lengthSqr()) * h);
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
            l = derivative.cross(r) / derivative.length();
        }
        return Vec2<Vec2Mode::SL, T>{s, l};
    }

    [[using gnu: pure]]
    Vec2<Vec2Mode::SL, T> inverse(Vec2<Vec2Mode::XY, T> point, T start_s, T end_s) const noexcept {
        if (start_s > end_s) {
            spdlog::warn(
                "Invalid argument issue detected! The start_s should always be less than end_s: "
                "start_s = {0:.6f} while end_s = {1:.6f}.",
                start_s, end_s
            );
            std::swap(start_s, end_s);
        }
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<Vec2Mode::XY, T>>& anchor_points = anchorPoints();
        const std::vector<Vec2<Vec2Mode::XY, T>>& ddys = vec2_of_s_.ddys();
        const std::vector<Vec2<Vec2Mode::XY, T>>& d4ys = vec2_of_s_.d4ys();
        const std::size_t istart =
            std::upper_bound(arc_lengths.cbegin(), arc_lengths.cend(), start_s) -
            arc_lengths.cbegin() - 1;
        const std::size_t iend = std::upper_bound(arc_lengths.cbegin(), arc_lengths.cend(), end_s) -
                                 arc_lengths.cbegin();
        const std::size_t pos =
            closetUpperElement(
                anchor_points.cbegin() + istart, anchor_points.cend() + iend, point
            ) -
            anchor_points.cbegin();
        Vec2<Vec2Mode::XY, T> derivative, derivative2, r;
        T h, derivative_length, s, l;
        if (pos == istart) {
            h = arc_lengths[1] - arc_lengths[0];
            r = point - anchor_points[0];
            derivative = quinerpd(
                anchor_points[0], anchor_points[1], ddys[0], ddys[1], d4ys[0], d4ys[1], 0.0, h
            );
            derivative_length = derivative.length();
            s = arc_lengths.front() + derivative.dot(r) / derivative_length;
            l = derivative.cross(r) / derivative_length;
        } else if (pos == iend) {
            h = arc_lengths[pos - 1] - arc_lengths[pos - 2];
            r = point - anchor_points[pos - 2];
            derivative = quinerpd(
                anchor_points[pos - 2], anchor_points[pos - 1], ddys[pos - 2], ddys[pos - 1],
                d4ys[pos - 2], d4ys[pos - 1], 1.0, h
            );
            derivative_length = derivative.length();
            s = arc_lengths.back() + derivative.dot(r) / derivative_length;
            l = derivative.cross(r) / derivative_length;
        } else {
            const Vec2<Vec2Mode::XY, T>& lower_anchor_point{anchor_points[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_anchor_point{anchor_points[pos]};
            const Vec2<Vec2Mode::XY, T>& lower_derivative2{ddys[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_derivative2{ddys[pos]};
            const Vec2<Vec2Mode::XY, T>& lower_derivative4{d4ys[pos - 1]};
            const Vec2<Vec2Mode::XY, T>& upper_derivative4{d4ys[pos]};
            T ratio = 0.5;
            h = vec2_of_s_.ts()[pos] - vec2_of_s_.ts()[pos - 1];
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
            for (std::size_t num_iter = 3; num_iter; --num_iter) {
                ratio -= r.dot(derivative) / ((r.dot(derivative2) - derivative.lengthSqr()) * h);
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
            l = derivative.cross(r) / derivative.length();
        }
        return Vec2<Vec2Mode::SL, T>{s, l};
    }

    [[using gnu: pure, always_inline]]
    T minS() const noexcept {
        return vec2_of_s_.minT();
    }

    [[using gnu: pure, always_inline]]
    T maxS() const noexcept {
        return vec2_of_s_.maxT();
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& arcLengths() const noexcept {
        return vec2_of_s_.ts();
    }

    [[using gnu: pure, always_inline]]
    const std::vector<Vec2<Vec2Mode::XY, T>>& anchorPoints() const noexcept {
        return vec2_of_s_.ys();
    }

  private:
    [[using gnu: const, always_inline]]
    static T calcArcLength(
        Vec2<Vec2Mode::XY, T> start, Vec2<Vec2Mode::XY, T> end, Vec2<Vec2Mode::XY, T> ddstart,
        Vec2<Vec2Mode::XY, T> ddend, Vec2<Vec2Mode::XY, T> d4start, Vec2<Vec2Mode::XY, T> d4end,
        T scale
    ) noexcept {
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
            0.95308992296933199639};
        constexpr std::array<T, 5> kGaussLegendreWeights{
            0.11846344252809454375, 0.23931433524968323402, 128.0 / 450.0, 0.23931433524968323402,
            0.11846344252809454375};
        T result{0.0};
        for (std::size_t i = 0; i < kNumGaussLegendreKnots; ++i) {
            result +=
                quinerpd(start, end, ddstart, ddend, d4start, d4end, kGaussLegendreKnots[i], scale)
                    .length() *
                kGaussLegendreWeights[i];
        }
        result *= scale;
        return result;
    }

    [[using gnu: always_inline]] explicit PiecewiseQuinticCurve2(
        PiecewiseQuinticFunction1<math::Vec2<Vec2Mode::XY, T>, T> vec2_of_s
    ) noexcept
        : vec2_of_s_(std::move(vec2_of_s)) {}

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& vec2_of_s_;
        return;
    }

    PiecewiseQuinticFunction1<Vec2<Vec2Mode::XY, T>, T> vec2_of_s_;
};

template <typename T>
struct isCurve<PiecewiseQuinticCurve2<T>> final {
    static constexpr bool value = true;
};

using PiecewiseQuinticCurve2f = PiecewiseQuinticCurve2<float>;
using PiecewiseQuinticCurve2d = PiecewiseQuinticCurve2<double>;

} // namespace math
} // namespace tiny_pnc
