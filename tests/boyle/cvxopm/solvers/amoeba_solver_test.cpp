/**
 * @file amoeba_solver_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-11-16
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/cvxopm/solvers/amoeba_solver.hpp"

#include <utility>

#include "boyle/cvxopm/problems/dense_problem.hpp"
#include "boyle/math/dense/matrixx.hpp"
#include "boyle/math/dense/vectorx.hpp"
#include "boyle/math/mdfunctions/quadratic_mdfunction.hpp"
#include "boyle/math/mdfunctions/rosenbrock_function.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::cvxopm {

TEST_CASE_TEMPLATE("QuadraticTest", T, float, double, long double) {
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

    const AmoebaSolver<T> amoeba_solver{.settings{.eps_abs{1E-6}}};

    ::boyle::math::pmr::VectorX<T> x0({0.5, 0.5});

    const auto [result, info] = amoeba_solver.solve(dense_problem, std::move(x0));

    CHECK_LE(info.iter, 23);
    CHECK_EQ(result.prim_vars[0], doctest::Approx(0.0).epsilon(9E-4));
    CHECK_EQ(result.prim_vars[1], doctest::Approx(2.0).epsilon(4E-4));
}

TEST_CASE_TEMPLATE("RosenBrockTest", T, float, double, long double) {
    const DenseProblem<T> dense_problem(
        ::boyle::math::makeMdFunctionProxy(
            ::boyle::math::RosenbrockFunction<::boyle::math::pmr::VectorX<T>>(1.0, 100.0)
        )
    );

    const AmoebaSolver<T> amoeba_solver{.settings{.eps_abs{1E-6}}};

    ::boyle::math::pmr::VectorX<T> x0({0.5, 0.5});

    const auto [result, info] = amoeba_solver.solve(dense_problem, std::move(x0));

    CHECK_LE(info.iter, 140);
    CHECK_EQ(result.prim_vars[0], doctest::Approx(1.0).epsilon(4E-2));
    CHECK_EQ(result.prim_vars[1], doctest::Approx(1.0).epsilon(7E-2));
}

} // namespace boyle::cvxopm
