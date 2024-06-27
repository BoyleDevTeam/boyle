/**
 * @file function1_proxy_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/functions/function1_proxy.hpp"

#include <vector>

#include "boyle/math/functions/piecewise_functions/piecewise_cubic_function1.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_linear_function1.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_quintic_function1.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Polymorphism") {
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
    PiecewiseLinearFunction1d linear_func{ts, ys};
    PiecewiseCubicFunction1d cubic_func{ts, ys};
    PiecewiseQuinticFunction1d quintic_func{ts, ys};

    std::vector<Function1Proxy<double, double>> function1s;
    function1s.emplace_back(makeFunction1Proxy(linear_func));
    function1s.emplace_back(makeFunction1Proxy(std::move(cubic_func)));
    function1s.emplace_back(makeFunction1Proxy(std::move(quintic_func)));

    for (const auto& func : function1s) {
        CHECK_EQ(func->minT(), doctest::Approx(-2.0).epsilon(kEpsilon));
        CHECK_EQ(func->maxT(), doctest::Approx(2.0).epsilon(kEpsilon));

        for (double t : ts) {
            CHECK_EQ(func->eval(t), exact_func(t));
        }

        for (double t : linspace(-2.0, 2.0, 31)) {
            CHECK_EQ(func->eval(t), doctest::Approx(exact_func(t)).epsilon(1E-2));
            CHECK_EQ(func->derivative(t), doctest::Approx(exact_dfunc(t)).epsilon(1E-1));
        }

        CHECK_EQ(
            func->integral(-2.0, 2.0), doctest::Approx(exact_integral(-2.0, 2.0)).epsilon(1E-2)
        );
        CHECK_EQ(
            func->integral(-1.46, -0.25),
            doctest::Approx(exact_integral(-1.46, -0.25)).epsilon(1E-3)
        );
        CHECK_EQ(
            func->integral(0.153, 1.97), doctest::Approx(exact_integral(0.153, 1.97)).epsilon(1E-3)
        );
    }
}

} // namespace boyle::math
