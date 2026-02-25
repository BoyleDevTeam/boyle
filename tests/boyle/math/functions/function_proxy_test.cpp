/**
 * @file function_proxy_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/functions/function_proxy.hpp"

#include <vector>

#include "boyle/math/functions/piecewise_cubic_function.hpp"
#include "boyle/math/functions/piecewise_linear_function.hpp"
#include "boyle/math/functions/piecewise_quintic_function.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Polymorphism") {
    constexpr auto exact_func = [](double t) noexcept -> double {
        return 0.45 + 5.3 * t - 1.3 * t * t + 0.65 * t * t * t + 0.075 * t * t * t * t -
               0.0014 * t * t * t * t * t;
    };
    constexpr auto exact_dfunc = [](double t) noexcept -> double {
        return 5.3 - 2.6 * t + 1.95 * t * t + 0.3 * t * t * t - 0.007 * t * t * t * t;
    };
    constexpr auto exact_integral = [](double lo, double up) noexcept -> double {
        return 0.45 * (up - lo) + 2.65 * (up * up - lo * lo) -
               1.3 / 3.0 * (up * up * up - lo * lo * lo) +
               0.1625 * (up * up * up * up - lo * lo * lo * lo) +
               0.015 * (up * up * up * up * up - lo * lo * lo * lo * lo) -
               0.0014 / 6.0 * (up * up * up * up * up * up - lo * lo * lo * lo * lo * lo);
    };

    constexpr std::size_t kNumPoints{41};
    constexpr double kStart{-2.0};
    constexpr double kEnd{2.0};

    const std::vector<double> ts = linspace(kStart, kEnd, kNumPoints);
    std::vector<double> ys;
    ys.reserve(kNumPoints);
    for (double t : ts) {
        ys.emplace_back(exact_func(t));
    }
    PiecewiseLinearFunction<double> linear_func{ts, ys};
    PiecewiseCubicFunction<double> cubic_func{ts, ys};
    PiecewiseQuinticFunction<double> quintic_func{ts, ys};

    std::vector<FunctionProxy<double>> function1s;
    function1s.emplace_back(makeFunctionProxy(linear_func));
    function1s.emplace_back(makeFunctionProxy(std::move(cubic_func)));
    function1s.emplace_back(makeFunctionProxy(std::move(quintic_func)));

    for (const auto& func : function1s) {
        CHECK_EQ(func->minT(), doctest::Approx(kStart).epsilon(kEpsilon));
        CHECK_EQ(func->maxT(), doctest::Approx(kEnd).epsilon(kEpsilon));

        for (double t : ts) {
            CHECK_EQ(func->eval(t), exact_func(t));
        }

        for (double t : linspace(kStart, kEnd, 31)) {
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
