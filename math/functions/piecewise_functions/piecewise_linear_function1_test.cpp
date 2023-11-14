/**
 * @file piecewise_linear_function_test1.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-11
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "math/functions/piecewise_functions/piecewise_linear_function1.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "math/utils.hpp"
#include "math/vec2.hpp"

namespace {

bool plot_graph;

} // namespace

namespace tiny_pnc {
namespace math {

TEST_CASE("TwoPointsLine") {
    const Vec2d a{0.0, 0.2};
    const Vec2d b{1.0, 0.4};
    auto exact_func = [](double t) -> double { return 0.2 * t + 0.2; };

    std::vector<double> ts = {a.x, b.x};
    std::vector<double> ys = {a.y, b.y};

    PiecewiseLinearFunction1d func{ts, ys};

    for (double t : linspace(-0.1, 1.1, 121)) {
        CHECK_EQ(func(t), doctest::Approx(exact_func(t)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(t), 0.2);
    }
}

TEST_CASE("QuinticPolynomial") {
    auto exact_func = [](double t) -> double {
        return 0.45 + 5.3 * t - 1.3 * t * t + 0.65 * t * t * t + 0.075 * t * t * t * t -
               0.0014 * t * t * t * t * t;
    };
    auto exact_dfunc = [](double t) -> double {
        return 5.3 - 2.6 * t + 1.95 * t * t + 0.3 * t * t * t - 0.007 * t * t * t * t;
    };
    auto exact_integral = [](double lo, double up) -> double {
        return 0.45 * (up - lo) + 2.65 * (up * up - lo * lo) -
               1.3 / 3.0 * (up * up * up - lo * lo * lo) +
               0.1625 * (up * up * up * up - lo * lo * lo * lo) +
               0.015 * (up * up * up * up * up - lo * lo * lo * lo * lo) -
               0.0014 / 6.0 * (up * up * up * up * up * up - lo * lo * lo * lo * lo * lo);
    };

    std::vector<double> ts = linspace(-2.0, 2.0, 41);
    std::vector<double> ys;
    ys.reserve(41);
    for (double t : ts) {
        ys.emplace_back(exact_func(t));
    }

    PiecewiseLinearFunction1d func{ts, ys};

    CHECK_EQ(func.minT(), doctest::Approx(-2.0).epsilon(kEpsilon));
    CHECK_EQ(func.maxT(), doctest::Approx(2.0).epsilon(kEpsilon));
    CHECK_EQ(func.minY(), doctest::Approx(exact_func(-2.0)).epsilon(kEpsilon));
    CHECK_EQ(func.maxY(), doctest::Approx(exact_func(2.0)).epsilon(kEpsilon));

    for (double t : ts) {
        CHECK_EQ(func(t), exact_func(t));
    }

    for (double t : linspace(-2.0, 2.0, 31)) {
        CHECK_EQ(func(t), doctest::Approx(exact_func(t)).epsilon(1E-2));
        CHECK_EQ(func.derivative(t), doctest::Approx(exact_dfunc(t)).epsilon(1E-1));
    }

    CHECK_EQ(func.integral(-2.0, 2.0), doctest::Approx(exact_integral(-2.0, 2.0)).epsilon(1E-2));
    CHECK_EQ(
        func.integral(-1.46, -0.25), doctest::Approx(exact_integral(-1.46, -0.25)).epsilon(1E-3)
    );
    CHECK_EQ(
        func.integral(0.153, 1.97), doctest::Approx(exact_integral(0.153, 1.97)).epsilon(1E-3)
    );

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = tiny_pnc::math::linspace(-4.0, 4.0, 801);
        {
            auto fig = figure();
            fig->title("func");
            hold(on);
            grid(on);
            plot(func.ts(), transform(func.ts(), exact_func), ".");
            plot(plot_ts, transform(plot_ts, func), "-");
            fig->show();
        }
        {
            auto fig = figure();
            fig->title("dfunc");
            hold(on);
            grid(on);
            plot(func.ts(), transform(func.ts(), exact_dfunc), ".");
            auto dfunc = [&func](double t) -> double { return func.derivative(t); };
            plot(plot_ts, transform(plot_ts, dfunc), "-");
            fig->show();
        }
    }
}

TEST_CASE("Serialization") {
    auto exact_func = [](double t) -> double {
        return 0.45 + 5.3 * t - 1.3 * t * t + 0.65 * t * t * t;
    };

    std::vector<double> ts = linspace(-2.0, 2.0, 41);
    std::vector<double> ys;
    for (double t : ts) {
        ys.emplace_back(exact_func(t));
    }

    PiecewiseLinearFunction1d func{ts, ys};

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << func;

    PiecewiseLinearFunction1d other_func;

    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_func;

    for (double t : ts) {
        CHECK_EQ(other_func(t), exact_func(t));
    }
}

} // namespace math
} // namespace tiny_pnc

int main(int argc, const char* argv[]) {
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
