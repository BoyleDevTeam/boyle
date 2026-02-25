/**
 * @file qp_problem_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-20
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/cvxopm/problems/qp_problem.hpp"

#include <sstream>
#include <vector>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "boyle/math/utils.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

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
class QpProblemTestFixture {
  public:
    QpProblemTestFixture(const QpProblemTestFixture& other) noexcept = delete;
    auto operator=(const QpProblemTestFixture& other) noexcept -> QpProblemTestFixture& = delete;
    QpProblemTestFixture(QpProblemTestFixture&& other) noexcept = delete;
    auto operator=(QpProblemTestFixture&& other) noexcept -> QpProblemTestFixture& = delete;
    ~QpProblemTestFixture() noexcept = default;

    QpProblemTestFixture() noexcept : qp_problem(2, 3) {
        qp_problem.addQuadCostTerm(0, 0, 2.0);
        qp_problem.addQuadCostTerm(1, 1, 1.0);
        qp_problem.addQuadCostTerm(0, 1, 1.0);
        qp_problem.addLinCostTerm(0, 1.0);
        qp_problem.addLinCostTerm(1, 1.0);
        qp_problem.updateConstrainTerm(0, {{0, 1.0}, {1, 1.0}}, 1.0, 1.0);
        qp_problem.updateConstrainTerm(1, {{0, 1.0}}, 0.0, 0.7);
        qp_problem.updateConstrainTerm(2, {{1, 1.0}}, 0.0, 0.7);
    }

  protected:
    QpProblem<double, int> qp_problem;
};

TEST_CASE_FIXTURE(QpProblemTestFixture, "Basic") {
    std::vector<double> state_vec;
    double exact_cost{0.0};
    SUBCASE("state_0") {
        state_vec = std::vector<double>{0.0, 0.0};
        exact_cost = 0.0;
        CHECK_FALSE(qp_problem.validate(state_vec));
    }

    SUBCASE("state_1") {
        state_vec = std::vector<double>{0.462, 0.538};
        exact_cost = state_vec[0] * state_vec[0] * 2.0 + state_vec[1] * state_vec[1] +
                     state_vec[0] * state_vec[1] + state_vec[0] + state_vec[1];
        CHECK(qp_problem.validate(state_vec));
    }

    SUBCASE("state_2") {
        state_vec = std::vector<double>{0.462, 0.538 + ::boyle::math::kEpsilon};
        exact_cost = state_vec[0] * state_vec[0] * 2.0 + state_vec[1] * state_vec[1] +
                     state_vec[0] * state_vec[1] + state_vec[0] + state_vec[1];
        CHECK_FALSE(qp_problem.validate(state_vec));
    }

    SUBCASE("state_3") {
        state_vec = std::vector<double>{0.462 - ::boyle::math::kEpsilon, 0.538};
        exact_cost = state_vec[0] * state_vec[0] * 2.0 + state_vec[1] * state_vec[1] +
                     state_vec[0] * state_vec[1] + state_vec[0] + state_vec[1];
        CHECK_FALSE(qp_problem.validate(state_vec));
    }

    CHECK_EQ(
        qp_problem.cost(state_vec), doctest::Approx(exact_cost).epsilon(::boyle::math::kEpsilon)
    );
}

TEST_CASE_FIXTURE(QpProblemTestFixture, "Serialization") {
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << qp_problem;

    QpProblem<double, int> other_problem;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_problem;

    const std::vector<double> state_vec{0.462, 0.538};
    const double exact_cost = state_vec[0] * state_vec[0] * 2.0 + state_vec[1] * state_vec[1] +
                              state_vec[0] * state_vec[1] + state_vec[0] + state_vec[1];
    CHECK(other_problem.validate(state_vec));
    CHECK_EQ(
        other_problem.cost(state_vec), doctest::Approx(exact_cost).epsilon(::boyle::math::kEpsilon)
    );
}

} // namespace boyle::cvxopm
