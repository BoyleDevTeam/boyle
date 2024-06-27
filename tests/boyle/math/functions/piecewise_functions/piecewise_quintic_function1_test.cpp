/**
 * @file piecewise_quintic_function_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-22
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/functions/piecewise_functions/piecewise_quintic_function1.hpp"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "boyle/math/utils.hpp"

namespace {

bool plot_graph{false};

auto createFigureHandle() noexcept -> matplot::figure_handle {
    matplot::figure_handle fig = matplot::figure();
    fig->size(1600, 1000);
    return fig;
}

auto createAxesHandles(matplot::figure_handle fig) noexcept -> std::vector<matplot::axes_handle> {
    std::vector<matplot::axes_handle> axes;
    axes.reserve(6);
    for (std::size_t i = 0; i < 6; ++i) {
        axes.emplace_back(fig->add_subplot(2, 3, i));
        axes[i]->grid(matplot::on);
        axes[i]->hold(matplot::on);
        axes[i]->font_size(12.0);
        axes[i]->x_axis().label_font_size(12.0);
        axes[i]->y_axis().label_font_size(12.0);
    }
    axes[0]->title("func");
    axes[1]->title("dfunc");
    axes[2]->title("ddfunc");
    axes[3]->title("d3func");
    axes[4]->title("d4func");
    axes[5]->title("d5func");
    return axes;
}

} // namespace

namespace boyle::math {

TEST_CASE("QuinticPolynomial") {
    auto exact_func = [](double t) -> double {
        return 0.45 + 5.3 * t - 1.3 * t * t + 0.65 * t * t * t + 0.075 * t * t * t * t -
               0.0014 * t * t * t * t * t;
    };
    auto exact_dfunc = [](double t) -> double {
        return 5.3 - 2.6 * t + 1.95 * t * t + 0.3 * t * t * t - 0.007 * t * t * t * t;
    };
    auto exact_ddfunc = [](double t) -> double {
        return -2.6 + 3.9 * t + 0.9 * t * t - 0.028 * t * t * t;
    };
    auto exact_d3func = [](double t) -> double { return 3.9 + 1.8 * t - 0.084 * t * t; };
    auto exact_d4func = [](double t) -> double { return 1.8 - 0.168 * t; };
    auto exact_d5func = []([[maybe_unused]] double t) -> double { return -0.168; };
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

    PiecewiseQuinticFunction1d::BoundaryMode b01{}, b02{}, bf1{}, bf2{};

    SUBCASE("Boundary1212") {
        b01.order = 1;
        b01.derivative = exact_dfunc(ts.front());
        b02.order = 2;
        b02.derivative = exact_ddfunc(ts.front());
        bf1.order = 1;
        bf1.derivative = exact_dfunc(ts.back());
        bf2.order = 2;
        bf2.derivative = exact_ddfunc(ts.back());
    }
    SUBCASE("Boundary2424") {
        b01.order = 2;
        b01.derivative = exact_ddfunc(ts.front());
        b02.order = 4;
        b02.derivative = exact_d4func(ts.front());
        bf1.order = 2;
        bf1.derivative = exact_ddfunc(ts.back());
        bf2.order = 4;
        bf2.derivative = exact_d4func(ts.back());
    }
    SUBCASE("Boundary1313") {
        b01.order = 1;
        b01.derivative = exact_dfunc(ts.front());
        b02.order = 3;
        b02.derivative = exact_d3func(ts.front());
        bf1.order = 1;
        bf1.derivative = exact_dfunc(ts.back());
        bf2.order = 3;
        bf2.derivative = exact_d3func(ts.back());
    }
    SUBCASE("Boundary2323") {
        b01.order = 2;
        b01.derivative = exact_ddfunc(ts.front());
        b02.order = 3;
        b02.derivative = exact_d3func(ts.front());
        bf1.order = 2;
        bf1.derivative = exact_ddfunc(ts.back());
        bf2.order = 3;
        bf2.derivative = exact_d3func(ts.back());
    }
    /*SUBCASE("Boundary3434") {
        b01.order = 3;
        b01.derivative = exact_d3func(ts.front());
        b02.order = 4;
        b02.derivative = exact_d4func(ts.front());
        bf1.order = 3;
        bf1.derivative = exact_d3func(ts.back());
        bf2.order = 4;
        bf2.derivative = exact_d4func(ts.back());
    }*/
    SUBCASE("Boundary1414") {
        b01.order = 1;
        b01.derivative = exact_dfunc(ts.front());
        b02.order = 4;
        b02.derivative = exact_d4func(ts.front());
        bf1.order = 1;
        bf1.derivative = exact_dfunc(ts.back());
        bf2.order = 4;
        bf2.derivative = exact_d4func(ts.back());
    }

    PiecewiseQuinticFunction1d func{ts, ys, {b01, b02}, {bf1, bf2}};

    CHECK_EQ(func.minT(), doctest::Approx(-2.0).epsilon(kEpsilon));
    CHECK_EQ(func.maxT(), doctest::Approx(2.0).epsilon(kEpsilon));
    CHECK_EQ(func.minY(), doctest::Approx(exact_func(-2.0)).epsilon(kEpsilon));
    CHECK_EQ(func.maxY(), doctest::Approx(exact_func(2.0)).epsilon(kEpsilon));

    for (auto it = ts.cbegin() + 1; it != ts.cend() - 1; ++it) {
        CHECK_EQ(func(*it), exact_func(*it));
        CHECK_EQ(func.derivative(*it), doctest::Approx(exact_dfunc(*it)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(*it, 2), doctest::Approx(exact_ddfunc(*it)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(*it, 3), doctest::Approx(exact_d3func(*it)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(*it, 4), doctest::Approx(exact_d4func(*it)).epsilon(kEpsilon));
    }

    for (double t : linspace(-2.0, 2.0, 31)) {
        CHECK_EQ(func(t), doctest::Approx(exact_func(t)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(t), doctest::Approx(exact_dfunc(t)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(t, 2), doctest::Approx(exact_ddfunc(t)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(t, 3), doctest::Approx(exact_d3func(t)).epsilon(kEpsilon));
        CHECK_EQ(func.derivative(t, 4), doctest::Approx(exact_d4func(t)).epsilon(kEpsilon));
    }

    CHECK_EQ(
        func.integral(-2.0, 2.0), doctest::Approx(exact_integral(-2.0, 2.0)).epsilon(kEpsilon)
    );
    CHECK_EQ(
        func.integral(-1.46, -0.25), doctest::Approx(exact_integral(-1.46, -0.25)).epsilon(kEpsilon)
    );
    CHECK_EQ(
        func.integral(0.153, 1.97), doctest::Approx(exact_integral(0.153, 1.97)).epsilon(kEpsilon)
    );

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(-4.0, 4.0, 801);

        figure_handle fig = createFigureHandle();
        std::vector<axes_handle> axes = createAxesHandles(fig);
        axes[0]->plot(func.ts(), transform(func.ts(), exact_func), ".");
        axes[0]->plot(plot_ts, transform(plot_ts, func), "-");
        axes[1]->plot(func.ts(), transform(func.ts(), exact_dfunc), ".");
        auto dfunc = [&func](double t) -> double { return func.derivative(t); };
        axes[1]->plot(plot_ts, transform(plot_ts, dfunc), "-");
        axes[2]->plot(func.ts(), transform(func.ts(), exact_ddfunc), ".");
        auto ddfunc = [&func](double t) -> double { return func.derivative(t, 2); };
        axes[2]->plot(plot_ts, transform(plot_ts, ddfunc), "-");
        axes[3]->plot(func.ts(), transform(func.ts(), exact_d3func), ".");
        auto d3func = [&func](double t) -> double { return func.derivative(t, 3); };
        axes[3]->plot(plot_ts, transform(plot_ts, d3func), "-");
        axes[4]->plot(func.ts(), transform(func.ts(), exact_d4func), ".");
        auto d4func = [&func](double t) -> double { return func.derivative(t, 4); };
        axes[4]->plot(plot_ts, transform(plot_ts, d4func), "-");
        axes[5]->plot(func.ts(), transform(func.ts(), exact_d5func), ".");
        auto d5func = [&func](double t) -> double { return func.derivative(t, 5); };
        axes[5]->plot(plot_ts, transform(plot_ts, d5func), "-");

        fig->show();
    }
}

TEST_CASE("Serialization") {
    auto exact_func = [](double t) -> double {
        return 0.45 + 5.3 * t - 1.3 * t * t + 0.65 * t * t * t + 0.075 * t * t * t * t -
               0.0014 * t * t * t * t * t;
    };

    std::vector<double> ts = linspace(-2.0, 2.0, 41);
    std::vector<double> ys;
    for (double t : ts) {
        ys.emplace_back(exact_func(t));
    }

    PiecewiseQuinticFunction1d func{ts, ys};

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << func;

    PiecewiseQuinticFunction1d other_func;

    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_func;

    for (double t : ts) {
        CHECK_EQ(other_func(t), exact_func(t));
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
