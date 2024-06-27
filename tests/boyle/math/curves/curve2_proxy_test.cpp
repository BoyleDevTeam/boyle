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

#include "boyle/math/curves/piecewise_curves/piecewise_cubic_curve2.hpp"
#include "boyle/math/curves/piecewise_curves/piecewise_linear_curve2.hpp"
#include "boyle/math/curves/piecewise_curves/piecewise_quintic_curve2.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Polymorphism") {

    auto exact_semi_circle = [](double theta) -> Vec2d {
        return Vec2d{2.0 * std::cos(theta), 2.0 * std::sin(theta)};
    };

    const std::size_t num_anchors = 101;
    const std::size_t num_intpls = 100;
    std::vector<double> anchor_thetas = linspace(0.0, M_PI, num_anchors);
    std::vector<double> intpl_thetas = linspace(0.015, M_PI - 0.015, num_intpls);

    std::vector<Vec2d> anchor_points;
    anchor_points.reserve(num_anchors);
    for (double theta : anchor_thetas) {
        anchor_points.push_back(exact_semi_circle(theta));
    }

    PiecewiseLinearCurve2d linear_curve{anchor_points};
    PiecewiseCubicCurve2d cubic_curve{anchor_points};
    PiecewiseQuinticCurve2d quintic_curve{anchor_points};

    std::vector<Curve2Proxy<double>> curves;
    curves.emplace_back(makeCurve2Proxy(linear_curve));
    curves.emplace_back(makeCurve2Proxy(std::move(cubic_curve)));
    curves.emplace_back(makeCurve2Proxy(std::move(quintic_curve)));

    for (const auto& curve : curves) {

        const std::vector<double>& arc_lengths{curve->arcLengths()};

        CHECK_EQ(curve->minS(), 0.0);
        CHECK_EQ(curve->maxS(), doctest::Approx(M_PI * 2.0).epsilon(1E-3));

        for (std::size_t i = 0; i < num_anchors; ++i) {
            const Vec2d p1 = curve->eval(arc_lengths[i]);
            const Vec2d p2 = exact_semi_circle(anchor_thetas[i]);
            CHECK_EQ(p1.x(), doctest::Approx(p2.x()).epsilon(kEpsilon));
            CHECK_EQ(p1.y(), doctest::Approx(p2.y()).epsilon(kEpsilon));
        }

        for (double theta : intpl_thetas) {
            const Vec2d point = exact_semi_circle(theta);
            const SlDupletd sl = curve->inverse(point);
            CHECK_EQ(sl.s(), doctest::Approx(2.0 * theta).epsilon(1E-3));
            CHECK_EQ(sl.l(), doctest::Approx(0.0).epsilon(1E-3));
        }
    }
}

} // namespace boyle::math
