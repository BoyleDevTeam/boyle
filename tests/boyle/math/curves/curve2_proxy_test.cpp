/**
 * @file curve2_proxy_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/curves/curve2_proxy.hpp"

#include <vector>

#include "boyle/math/curves/piecewise_cubic_curve.hpp"
#include "boyle/math/curves/piecewise_linear_curve.hpp"
#include "boyle/math/curves/piecewise_quintic_curve.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Polymorphism") {
    const auto exact_semi_circle = [](double theta) noexcept -> Vec2d {
        return Vec2d{2.0 * std::cos(theta), 2.0 * std::sin(theta)};
    };

    constexpr std::size_t kNumAnchors = 101;
    constexpr double kStart = 0.0;
    constexpr double kEnd = M_PI;
    constexpr double kStep = (kEnd - kStart) / (kNumAnchors - 1);

    std::vector<Vec2d> anchor_points(kNumAnchors);
    std::ranges::generate(
        anchor_points,
        [t = kStart - kStep, h = kStep]() mutable noexcept -> Vec2d {
            t += h;
            return Vec2d{2.0 * std::cos(t), 2.0 * std::sin(t)};
        }
    );

    PiecewiseLinearCurve2d linear_curve{anchor_points};
    PiecewiseCubicCurve2d cubic_curve{anchor_points};
    PiecewiseQuinticCurve2d quintic_curve{anchor_points};

    std::vector<Curve2Proxy<Vec2d>> curves;
    curves.emplace_back(makeCurve2Proxy(linear_curve));
    curves.emplace_back(makeCurve2Proxy(std::move(cubic_curve)));
    curves.emplace_back(makeCurve2Proxy(std::move(quintic_curve)));

    for (const auto& curve : curves) {
        const std::vector<double>& arc_lengths{curve->arcLengths()};

        CHECK_EQ(curve->minS(), 0.0);
        CHECK_EQ(curve->maxS(), doctest::Approx(M_PI * 2.0).epsilon(1E-3));

        for (std::size_t i = 0; i < kNumAnchors; ++i) {
            const double theta{kStart + kStep * i};
            const Vec2d p1 = curve->eval(arc_lengths[i]);
            const Vec2d p2 = exact_semi_circle(theta);
            CHECK_EQ(p1.x, doctest::Approx(p2.x).epsilon(kEpsilon));
            CHECK_EQ(p1.y, doctest::Approx(p2.y).epsilon(kEpsilon));
        }

        for (std::size_t i = 1; i < kNumAnchors; ++i) {
            const double theta{kStart + kStep * (i - 0.5)};
            const Vec2d point = exact_semi_circle(theta);
            const SlDupletd sl = curve->inverse(point);
            CHECK_EQ(sl.s, doctest::Approx(2.0 * theta).epsilon(1E-3));
            CHECK_EQ(sl.l, doctest::Approx(0.0).epsilon(1E-3));
        }
    }
}

} // namespace boyle::math
