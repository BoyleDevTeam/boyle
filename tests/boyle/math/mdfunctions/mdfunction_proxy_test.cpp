/**
 * @file mdfunction_proxy_test.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-14
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/math/mdfunctions/mdfunction_proxy.hpp"

#include "boyle/math/dense/detail/dense_generate_trait.hpp"
#include "boyle/math/dense/matrix.hpp"
#include "boyle/math/dense/matrixx.hpp"
#include "boyle/math/dense/vector.hpp"
#include "boyle/math/dense/vectorx.hpp"
#include "boyle/math/mdfunctions/linear_mdfunction.hpp"
#include "boyle/math/mdfunctions/quadratic_mdfunction.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE_TEMPLATE(
    "Polymorphism", T, Vector<float, 5>, Vector<double, 5>, VectorX<float>, VectorX<double>
) {
    T linear_coeffs(5);
    linear_coeffs[0] = 674.0;
    linear_coeffs[1] = -194.45;
    linear_coeffs[2] = 22.55;
    linear_coeffs[3] = 1.056;
    linear_coeffs[4] = -0.47656;

    ::boyle::math::detail::DenseGenerateTraitT<T> quadratic_coeffs(5, 5);
    quadratic_coeffs[0, 0] = 498.1031;
    quadratic_coeffs[0, 1] = -463.8805;
    quadratic_coeffs[0, 2] = 576.4003;
    quadratic_coeffs[0, 3] = -433.1227;
    quadratic_coeffs[0, 4] = -241.3050;
    quadratic_coeffs[1, 0] = -463.8805;
    quadratic_coeffs[1, 1] = -390.2700;
    quadratic_coeffs[1, 2] = -192.6189;
    quadratic_coeffs[1, 3] = 356.3896;
    quadratic_coeffs[1, 4] = -951.7053;
    quadratic_coeffs[2, 0] = 576.4003;
    quadratic_coeffs[2, 1] = -192.6189;
    quadratic_coeffs[2, 2] = -95.0645;
    quadratic_coeffs[2, 3] = -305.2551;
    quadratic_coeffs[2, 4] = -186.8198;
    quadratic_coeffs[3, 0] = -433.1227;
    quadratic_coeffs[3, 1] = 356.3896;
    quadratic_coeffs[3, 2] = -305.2551;
    quadratic_coeffs[3, 3] = -948.3788;
    quadratic_coeffs[3, 4] = 69.7871;
    quadratic_coeffs[4, 0] = -241.3050;
    quadratic_coeffs[4, 1] = -951.7053;
    quadratic_coeffs[4, 2] = -186.8198;
    quadratic_coeffs[4, 3] = 69.7871;
    quadratic_coeffs[4, 4] = -919.1390;

    LinearMdFunction<T> linear_mdfunction{1.264, linear_coeffs};
    QuadraticMdFunction<T> quadratic_mdfunction(1.264, linear_coeffs, quadratic_coeffs);

    std::vector<MdFunctionProxy<T>> mdfunctions{};
    mdfunctions.emplace_back(makeMdFunctionProxy(linear_mdfunction));
    mdfunctions.emplace_back(makeMdFunctionProxy(std::move(quadratic_mdfunction)));

    T x(5);
    x[0] = 365.0;
    x[1] = -586.44;
    x[2] = 0.0;
    x[3] = 35.46;
    x[4] = -75.468;

    CHECK_EQ(mdfunctions[0]->num_dimensions(), 5);
    CHECK_EQ(mdfunctions[1]->num_dimensions(), 5);

    CHECK_EQ(mdfunctions[0]->eval(x), doctest::Approx(x.dot(linear_coeffs) + 1.264).epsilon(1E-6));
    CHECK_EQ(
        mdfunctions[1]->eval(x),
        doctest::Approx(x.dot(quadratic_coeffs).dot(x) * 0.5 + x.dot(linear_coeffs) + 1.264)
            .epsilon(1E-6)
    );

    CHECK(mdfunctions[0]->gradient(x).identicalTo(linear_coeffs, 1E-6));
    CHECK(mdfunctions[1]->gradient(x).identicalTo(linear_coeffs + x.dot(quadratic_coeffs), 1E-6));

    CHECK_FALSE(mdfunctions[0]->has_extrema(x));
    CHECK_FALSE(mdfunctions[1]->has_extrema(x));
}

} // namespace boyle::math
