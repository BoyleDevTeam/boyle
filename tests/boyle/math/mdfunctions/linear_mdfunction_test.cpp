/**
 * @file linear_mdfunction_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-09-07
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/math/mdfunctions/linear_mdfunction.hpp"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "boyle/math/dense/vector.hpp"
#include "boyle/math/dense/vectorx.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE_TEMPLATE(
    "Basic", T, Vector<float, 5>, Vector<double, 5>, VectorX<float>, VectorX<double>
) {
    T linear_coeffs(5);
    linear_coeffs[0] = 674.0;
    linear_coeffs[1] = -194.45;
    linear_coeffs[2] = 22.55;
    linear_coeffs[3] = 1.056;
    linear_coeffs[4] = -0.47656;

    LinearMdFunction<T> linear_mdfunction{1.264, linear_coeffs};

    CHECK_EQ(linear_mdfunction.num_dimensions(), 5);

    CHECK_EQ(linear_mdfunction.bias(), doctest::Approx(1.264).epsilon(1E-6));
    for (std::size_t i{0}; i < 5; ++i) {
        CHECK_EQ(linear_mdfunction.linearCoeff(i), doctest::Approx(linear_coeffs[i]).epsilon(1E-6));
    }

    linear_mdfunction.updateBias(103.984);
    linear_mdfunction.updateLinearCoeff(3, 2.475);

    CHECK_NE(linear_mdfunction.bias(), doctest::Approx(1.264).epsilon(1E-6));
    CHECK_NE(linear_mdfunction.linearCoeff(3), doctest::Approx(linear_coeffs[3]).epsilon(1E-6));

    CHECK_EQ(linear_mdfunction.bias(), doctest::Approx(103.984).epsilon(1E-6));
    CHECK_EQ(linear_mdfunction.linearCoeff(3), doctest::Approx(2.475).epsilon(1E-6));

    linear_mdfunction.updateBias(1.264);
    linear_mdfunction.updateLinearCoeff(3, 1.056);

    T x(5);
    x[0] = 365.0;
    x[1] = -586.44;
    x[2] = 0.0;
    x[3] = 35.46;
    x[4] = -75.468;

    CHECK_EQ(linear_mdfunction(x), doctest::Approx(linear_coeffs.dot(x) + 1.264).epsilon(1E-6));

    const T gradient{linear_mdfunction.gradient(x)};
    for (std::size_t i{0}; i < 5; ++i) {
        CHECK_EQ(gradient[i], doctest::Approx(linear_coeffs[i]).epsilon(1E-6));
    }

    for (std::size_t i{0}; i < 5; ++i) {
        CHECK_EQ(linear_mdfunction.gradient(x, i), doctest::Approx(linear_coeffs[i]).epsilon(1E-6));
    }
}

TEST_CASE_TEMPLATE(
    "Serialization", T, Vector<float, 5>, Vector<double, 5>, VectorX<float>, VectorX<double>
) {
    T linear_coeffs(5);
    linear_coeffs[0] = 674.0;
    linear_coeffs[1] = -194.45;
    linear_coeffs[2] = 22.55;
    linear_coeffs[3] = 1.056;
    linear_coeffs[4] = -0.47656;

    LinearMdFunction<T> linear_mdfunction{1.264, linear_coeffs};

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << linear_mdfunction;

    LinearMdFunction<T> other_Linear_mdfunction;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_Linear_mdfunction;

    CHECK_EQ(linear_mdfunction.bias(), other_Linear_mdfunction.bias());
    for (std::size_t i{0}; i < 5; ++i) {
        CHECK_EQ(linear_mdfunction.linearCoeff(i), other_Linear_mdfunction.linearCoeff(i));
    }
}

} // namespace boyle::math
