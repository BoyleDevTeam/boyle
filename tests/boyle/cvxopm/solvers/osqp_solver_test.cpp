/**
 * @file osqp_solver_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-21
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/cvxopm/solvers/osqp_solver.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "boyle/common/utils/macros.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::cvxopm {

/**
 * @brief Here we construct a simple Quadratic programming model to test.
 *        cost:
 *            f(x0, x1) = 2 * x0^2 + x1^2 + x0*x1 + x0 + x1
 *        constraints:
 *            1 <= x0 + x1 <= 1
 *            0 <= x0 <= 0.7
 *            0 <= x1 <= 0.7
 *
 */
class OsqpSolverTestFixture {
  public:
    OsqpSolverTestFixture() : qp_problem(2, 3), solver() {
        solver.settings.polishing = 1;

        qp_problem.updateQuadCostTerm(0, 0, 2.0);
        qp_problem.updateQuadCostTerm(1, 1, 1.0);
        qp_problem.updateQuadCostTerm(0, 1, 1.0);
        qp_problem.updateLinCostTerm(0, 1.0);
        qp_problem.updateLinCostTerm(1, 1.0);
        qp_problem.updateConstrainTerm(0, {{0, 1.0}, {1, 1.0}}, 1.0, 1.0);
        qp_problem.updateConstrainTerm(1, {{0, 1.0}}, 0.0, 0.7);
        qp_problem.updateConstrainTerm(2, {{1, 1.0}}, 0.0, 0.7);

        // Run cold start
        std::pair<OsqpSolver::Result, OsqpSolver::Info> pair = solver.solve(qp_problem);
        result = std::move(pair.first);
        info = pair.second;

        CHECK_EQ(result.prim_vars[0], doctest::Approx(0.3).epsilon(1e-7));
        CHECK_EQ(result.prim_vars[1], doctest::Approx(0.7).epsilon(1e-7));
        CHECK_EQ(result.dual_vars[0], doctest::Approx(-2.9).epsilon(1e-7));
        CHECK_EQ(result.dual_vars[1], 0.0);
        CHECK_EQ(result.dual_vars[2], doctest::Approx(0.2).epsilon(1e-6));
        CHECK_EQ(info.iter, 25);
    }
    DISABLE_COPY_AND_MOVE(OsqpSolverTestFixture);
    ~OsqpSolverTestFixture() noexcept = default;

  protected:
    QpProblem<double, int> qp_problem{};
    OsqpSolver::Result result{};
    OsqpSolver::Info info{};
    OsqpSolver solver{};
};

TEST_CASE_FIXTURE(OsqpSolverTestFixture, "WarmStart") {
    const std::vector<double> prim_vars_0{0.3, 0.7};
    const std::vector<double> dual_vars_0{-2.9, 0.0, 0.2};

    const auto [result, info] = solver.solve(qp_problem, prim_vars_0, dual_vars_0);

    CHECK_EQ(result.prim_vars[0], doctest::Approx(0.3).epsilon(1e-7));
    CHECK_EQ(result.prim_vars[1], doctest::Approx(0.7).epsilon(1e-7));
    CHECK_EQ(result.dual_vars[0], doctest::Approx(-2.9).epsilon(1e-7));
    CHECK_EQ(result.dual_vars[1], 0.0);
    CHECK_EQ(result.dual_vars[2], doctest::Approx(0.2).epsilon(1e-6));
    CHECK_EQ(info.iter, 25);
}

TEST_SUITE("Serialization") {
    TEST_CASE_FIXTURE(OsqpSolverTestFixture, "OsqpSolver::Result") {
        std::ostringstream oss;
        boost::archive::binary_oarchive oa(oss);
        oa << result;

        OsqpSolver::Result other_result;
        std::istringstream iss(oss.str());
        boost::archive::binary_iarchive ia(iss);
        ia >> other_result;

        CHECK_EQ(other_result.prim_vars[0], doctest::Approx(0.3).epsilon(1e-7));
        CHECK_EQ(other_result.prim_vars[1], doctest::Approx(0.7).epsilon(1e-7));
        CHECK_EQ(other_result.dual_vars[0], doctest::Approx(-2.9).epsilon(1e-7));
        CHECK_EQ(other_result.dual_vars[1], 0.0);
        CHECK_EQ(other_result.dual_vars[2], doctest::Approx(0.2).epsilon(1e-6));
    }

    TEST_CASE_FIXTURE(OsqpSolverTestFixture, "OsqpSolver") {
        std::ostringstream oss;
        boost::archive::binary_oarchive oa(oss);
        oa << solver;

        OsqpSolver other_solver;
        std::istringstream iss(oss.str());
        boost::archive::binary_iarchive ia(iss);
        ia >> other_solver;

        CHECK_EQ(other_solver.settings.scaling, solver.settings.scaling);
        CHECK_EQ(other_solver.settings.warm_starting, solver.settings.warm_starting);
        CHECK_EQ(other_solver.settings.polishing, solver.settings.polishing);
        CHECK_EQ(other_solver.settings.rho, solver.settings.rho);
        CHECK_EQ(other_solver.settings.alpha, solver.settings.alpha);
        CHECK_EQ(other_solver.settings.sigma, solver.settings.sigma);
    }
}

} // namespace boyle::cvxopm
