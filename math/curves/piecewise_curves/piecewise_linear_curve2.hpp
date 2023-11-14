/**
 * @file piecewise_linear_curve2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-28
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <format>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "boost/serialization/access.hpp"
#include "spdlog/spdlog.h"

#include "common/utils/macros.hpp"
#include "math/curves/curve2.hpp"
#include "math/functions/piecewise_functions/piecewise_linear_function1.hpp"
#include "math/utils.hpp"
#include "math/vec2.hpp"

namespace tiny_pnc {
namespace math {

template <typename T>
class [[nodiscard]] PiecewiseLinearCurve2 final : public Curve2<PiecewiseLinearCurve2<T>, T> {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    static constexpr T kDuplicateCriterion{kEpsilon};

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseLinearCurve2);

    [[using gnu: flatten, leaf]] explicit PiecewiseLinearCurve2(
        std::vector<Vec2<Vec2Mode::XY, T>> anchor_points, T s0 = 0.0
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
        std::vector<T> arc_lengths;
        arc_lengths.reserve(size);
        arc_lengths.push_back(s0);
        for (typename std::vector<Vec2<Vec2Mode::XY, T>>::const_iterator it =
                 anchor_points.cbegin() + 1;
             it != anchor_points.cend(); ++it) {
            arc_lengths.push_back(arc_lengths.back() + (it - 1)->distanceTo(*it));
        }
        vec2_of_s_ = PiecewiseLinearFunction1<Vec2<Vec2Mode::XY, T>, T>{
            std::move(arc_lengths), std::move(anchor_points)};
    }

    ~PiecewiseLinearCurve2() noexcept override = default;

    [[using gnu: pure, always_inline]]
    Vec2<Vec2Mode::XY, T>
    operator()(T s) const noexcept {
        return vec2_of_s_(s);
    }

    [[using gnu: pure, always_inline]]
    Vec2<Vec2Mode::XY, T> tangent(T s) const noexcept {
        return normalize(vec2_of_s_.derivative(s));
    }

    [[using gnu: pure, flatten, leaf]]
    Vec2<Vec2Mode::SL, T> inverse(Vec2<Vec2Mode::XY, T> point) const noexcept {
        const std::vector<T>& arc_lengths = arcLengths();
        const std::vector<Vec2<Vec2Mode::XY, T>>& anchor_points = anchorPoints();
        const std::size_t pos =
            closetUpperElement(anchor_points.cbegin(), anchor_points.cend(), point) -
            anchor_points.cbegin();
        Vec2<Vec2Mode::XY, T> diff, r;
        T diff_length, s, l;
        if (pos == 0) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_length = diff.length();
            s = arc_lengths.front() + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
        } else if (pos == anchor_points.size()) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_length = diff.length();
            s = arc_lengths.back() + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_length = diff.length();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
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
        Vec2<Vec2Mode::XY, T> diff, r;
        T diff_length, s, l;
        if (pos == istart) {
            r = point - anchor_points[0];
            diff = anchor_points[1] - anchor_points[0];
            diff_length = diff.length();
            s = arc_lengths.front() + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
        } else if (pos == iend) {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos - 1] - anchor_points[pos - 2];
            diff_length = diff.length();
            s = arc_lengths.back() + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
        } else {
            r = point - anchor_points[pos - 1];
            diff = anchor_points[pos] - anchor_points[pos - 1];
            diff_length = diff.length();
            s = arc_lengths[pos - 1] + diff.dot(r) / diff_length;
            l = diff.cross(r) / diff_length;
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
    [[using gnu: always_inline]] explicit PiecewiseLinearCurve2(
        PiecewiseLinearFunction1<math::Vec2<Vec2Mode::XY, T>, T> vec2_of_s
    )
        : vec2_of_s_(std::move(vec2_of_s)) {}

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& vec2_of_s_;
        return;
    }

    PiecewiseLinearFunction1<Vec2<Vec2Mode::XY, T>, T> vec2_of_s_;
};

template <typename T>
struct isCurve<PiecewiseLinearCurve2<T>> final {
    static constexpr bool value = true;
};

using PiecewiseLinearCurve2f = PiecewiseLinearCurve2<float>;
using PiecewiseLinearCurve2d = PiecewiseLinearCurve2<double>;

} // namespace math
} // namespace tiny_pnc
