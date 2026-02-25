/**
 * @file dense_problem_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-11-15
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/cvxopm/problems/dense_problem.hpp"

#include <utility>

#include "boyle/math/dense/matrixx.hpp"
#include "boyle/math/dense/vectorx.hpp"
#include "boyle/math/mdfunctions/quadratic_mdfunction.hpp"
#include "boyle/math/mdfunctions/rosenbrock_function.hpp"
#include "boyle/math/utils.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::cvxopm {

TEST_CASE_TEMPLATE("QuadraticTest", T, float, double, long double) {
    auto exact_quadratic_mdfunction = [](const ::boyle::math::pmr::VectorX<T>& x) noexcept -> T {
        return x[0] * x[0] + x[0] * x[1] + x[1] * x[1] - x[0] * static_cast<T>(2.0) -
               x[1] * static_cast<T>(4.0) + static_cast<T>(1.0);
    };
    ::boyle::math::QuadraticMdFunction<::boyle::math::pmr::VectorX<T>> quadratic_mdfunction(2);

    quadratic_mdfunction.updateBias(1.0);
    quadratic_mdfunction.updateLinearCoeff(0, -2.0);
    quadratic_mdfunction.updateLinearCoeff(1, -4.0);
    quadratic_mdfunction.updateQuadraticCoeff(0, 0, 1.0);
    quadratic_mdfunction.updateQuadraticCoeff(0, 1, 1.0);
    quadratic_mdfunction.updateQuadraticCoeff(1, 1, 1.0);

    const DenseProblem<T> dense_problem(
        ::boyle::math::makeMdFunctionProxy(std::move(quadratic_mdfunction))
    );

    const std::vector<double> x0s = ::boyle::math::linspace(-2.0, 4.0, 400);
    const std::vector<double> x1s = ::boyle::math::linspace(-2.0, 4.0, 400);

    CHECK_EQ(dense_problem.num_variables(), 2);

    for (T x0 : x0s) {
        for (T x1 : x1s) {
            const ::boyle::math::pmr::VectorX<T> x{{x0, x1}};
            CHECK_EQ(
                dense_problem.cost(x), doctest::Approx(exact_quadratic_mdfunction(x)).epsilon(1E-5)
            );
        }
    }
}

TEST_CASE_TEMPLATE("RosenBrockTest", T, float, double, long double) {
    auto exact_rosenbrock_function = [](const ::boyle::math::pmr::VectorX<T>& x) noexcept -> T {
        return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };
    const DenseProblem<T> dense_problem(
        ::boyle::math::makeMdFunctionProxy(
            ::boyle::math::RosenbrockFunction<::boyle::math::pmr::VectorX<T>>(1.0, 100.0)
        )
    );

    const std::vector<double> x0s = ::boyle::math::linspace(-2.0, 2.0, 400);
    const std::vector<double> x1s = ::boyle::math::linspace(-1.0, 3.0, 400);

    CHECK_EQ(dense_problem.num_variables(), 2);

    for (T x0 : x0s) {
        for (T x1 : x1s) {
            const ::boyle::math::pmr::VectorX<T> x{{x0, x1}};
            CHECK_EQ(
                dense_problem.cost(x), doctest::Approx(exact_rosenbrock_function(x)).epsilon(1E-5)
            );
        }
    }
}

} // namespace boyle::cvxopm
