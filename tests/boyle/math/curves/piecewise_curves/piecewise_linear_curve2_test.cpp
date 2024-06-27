/**
 * @file piecewise_linear_curve2_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-29
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/curves/piecewise_curves/piecewise_linear_curve2.hpp"

#include <cmath>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace {

bool plot_graph{false};

} // namespace

namespace boyle::math {

TEST_CASE("SemiCircle") {
    auto exact_semi_circle = [](double theta) -> Vec2d {
        return Vec2d{2.0 * std::cos(theta), 2.0 * std::sin(theta)};
    };

    const std::size_t num_anchors = 1001;
    const std::size_t num_intpls = 1000;
    std::vector<double> anchor_thetas = linspace(0.0, M_PI, num_anchors);
    std::vector<double> intpl_thetas = linspace(0.015, M_PI - 0.015, num_intpls);

    std::vector<Vec2d> anchor_points;
    anchor_points.reserve(num_anchors);
    for (double theta : anchor_thetas) {
        anchor_points.push_back(exact_semi_circle(theta));
    }

    PiecewiseLinearCurve2d semi_circle{anchor_points};

    const std::vector<double>& arc_lengths{semi_circle.arcLengths()};

    CHECK_EQ(semi_circle.minS(), 0.0);
    CHECK_EQ(semi_circle.maxS(), doctest::Approx(M_PI * 2.0).epsilon(1E-3));

    for (std::size_t i = 0; i < num_anchors; ++i) {
        const Vec2d p1 = semi_circle(arc_lengths[i]);
        const Vec2d p2 = exact_semi_circle(anchor_thetas[i]);
        CHECK_EQ(p1.x(), doctest::Approx(p2.x()).epsilon(kEpsilon));
        CHECK_EQ(p1.y(), doctest::Approx(p2.y()).epsilon(kEpsilon));
    }

    for (double theta : intpl_thetas) {
        const Vec2d point = exact_semi_circle(theta);
        const SlDupletd sl = semi_circle.inverse(point);
        CHECK_EQ(sl.s(), doctest::Approx(2.0 * theta).epsilon(3E-4));
        CHECK_EQ(sl.l(), doctest::Approx(0.0).epsilon(3E-4));
    }

    const std::array<Vec2d, 5> test_points{
        Vec2d{-1.3, 1.8}, Vec2d{-0.2, 1.8}, Vec2d{1000.0, 0}, Vec2d{0.8, 1.9}, Vec2d{1.6, 0.3}
    };
    for (Vec2d point : test_points) {
        CHECK_EQ(2.0 - point.norm(), doctest::Approx(semi_circle.inverse(point).l()).epsilon(3E-4));
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> exact_xs, exact_ys, intpl_xs, intpl_ys;
        exact_xs.reserve(num_anchors);
        exact_ys.reserve(num_anchors);
        intpl_xs.reserve(num_intpls);
        intpl_ys.reserve(num_intpls);
        for (Vec2d anchor_point : anchor_points) {
            exact_xs.push_back(anchor_point.x());
            exact_ys.push_back(anchor_point.y());
        }
        for (double s : linspace(semi_circle.minS(), semi_circle.maxS(), num_intpls)) {
            const Vec2d intpl_point = semi_circle(s);
            intpl_xs.push_back(intpl_point.x());
            intpl_ys.push_back(intpl_point.y());
        }
        {
            figure_handle fig = figure();
            fig->size(1600, 900);
            fig->title("func");
            hold(on);
            grid(on);
            axis(on);
            plot(exact_xs, exact_ys, ".");
            plot(intpl_xs, intpl_ys, "-");
            fig->show();
        }
    }
}

TEST_CASE("LogarithmicSpiral") {
    auto exact_log_spiral = [](double theta) -> Vec2d {
        return Vec2d{
            0.5 * std::exp(0.2 * theta) * std::cos(theta),
            0.5 * std::exp(0.2 * theta) * std::sin(theta)
        };
    };

    const std::size_t num_anchors = 1001;
    const std::size_t num_intpls = 1000;
    std::vector<double> anchor_thetas = linspace(0.0, M_PI * 8.0, num_anchors);
    std::vector<double> intpl_thetas = linspace(0.01, M_PI * 8.0 - 0.01, num_intpls);

    std::vector<Vec2d> anchor_points;
    anchor_points.reserve(num_anchors);
    for (double theta : anchor_thetas) {
        anchor_points.push_back(exact_log_spiral(theta));
    }

    PiecewiseLinearCurve2d log_spiral{anchor_points};

    const std::vector<double>& arc_lengths{log_spiral.arcLengths()};

    for (std::size_t i = 0; i < num_anchors; ++i) {
        const Vec2d p1 = log_spiral(arc_lengths[i]);
        const Vec2d p2 = exact_log_spiral(anchor_thetas[i]);
        CHECK_EQ(p1.x(), doctest::Approx(p2.x()).epsilon(kEpsilon));
        CHECK_EQ(p1.y(), doctest::Approx(p2.y()).epsilon(kEpsilon));
    }

    for (double theta : intpl_thetas) {
        const Vec2d point = exact_log_spiral(theta);
        const SlDupletd sl = log_spiral.inverse(point);
        CHECK_EQ(sl.l(), doctest::Approx(0.0).epsilon(6E-3));
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> exact_xs, exact_ys, intpl_xs, intpl_ys;
        exact_xs.reserve(num_anchors);
        exact_ys.reserve(num_anchors);
        intpl_xs.reserve(10001);
        intpl_ys.reserve(10001);
        for (Vec2d anchor_point : anchor_points) {
            exact_xs.push_back(anchor_point.x());
            exact_ys.push_back(anchor_point.y());
        }
        for (double s : linspace(log_spiral.minS(), log_spiral.maxS(), 10001)) {
            const Vec2d intpl_point = log_spiral(s);
            intpl_xs.push_back(intpl_point.x());
            intpl_ys.push_back(intpl_point.y());
        }

        figure_handle fig = figure();
        fig->size(1600, 900);
        fig->title("func");
        hold(on);
        grid(on);
        axis(on);
        plot(exact_xs, exact_ys, ".");
        plot(intpl_xs, intpl_ys, "-");
        fig->show();
    }
}

TEST_CASE("Serialization") {
    auto exact_log_spiral = [](double theta) -> Vec2d {
        return Vec2d{
            0.5 * std::exp(0.2 * theta) * std::cos(theta),
            0.5 * std::exp(0.2 * theta) * std::sin(theta)
        };
    };

    const std::size_t num_anchors = 1001;
    std::vector<double> anchor_thetas = linspace(0.0, M_PI * 8.0, num_anchors);

    std::vector<Vec2d> anchor_points;
    anchor_points.reserve(num_anchors);
    for (double theta : anchor_thetas) {
        anchor_points.push_back(exact_log_spiral(theta));
    }

    PiecewiseLinearCurve2d log_spiral{anchor_points};

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << log_spiral;

    PiecewiseLinearCurve2d other_log_spiral;

    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_log_spiral;

    const std::vector<double>& arc_length = log_spiral.arcLengths();
    const std::vector<double>& other_arc_length = other_log_spiral.arcLengths();
    const std::vector<Vec2d>& other_anchor_points = other_log_spiral.anchorPoints();

    for (std::size_t i = 0; i < num_anchors; ++i) {
        CHECK_EQ(other_arc_length[i], doctest::Approx(arc_length[i]).epsilon(kEpsilon));
        CHECK_EQ(
            other_anchor_points[i].x(), doctest::Approx(anchor_points[i].x()).epsilon(kEpsilon)
        );
        CHECK_EQ(
            other_anchor_points[i].y(), doctest::Approx(anchor_points[i].y()).epsilon(kEpsilon)
        );
    }
}

} // namespace boyle::math

auto main(int argc, const char* argv[]) -> int {
    cxxopts::Options options(
        "piecewise_linear_function_test", "unit test of PiecewiseLinearFunction class"
    );
    options.add_options()(
        "plot-graph", "plot test graph", cxxopts::value<bool>()->default_value("false")
    );
    cxxopts::ParseResult result = options.parse(argc, argv);
    plot_graph = result["plot-graph"].as<bool>();
    doctest::Context context(argc, argv);
    return context.run();
}
