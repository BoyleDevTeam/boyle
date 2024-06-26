/**
 * @file cubic_interpolation_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/cubic_interpolation.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "boyle/math/utils.hpp"

namespace boyle::math {

TEST_CASE("CubicPolynomial") {
    auto func = [](double t) -> double {
        return 28.483 - 2.34 * t + 1.54 * t * t + 0.365 * t * t * t;
    };
    auto dfunc = [](double t) -> double { return -2.34 + 3.08 * t + 1.095 * t * t; };
    auto ddfunc = [](double t) -> double { return 3.08 + 2.19 * t; };

    const double lo = -23.576;
    const double up = 97.273;
    const double start = func(lo);
    const double end = func(up);
    const double ddstart = ddfunc(lo);
    const double ddend = ddfunc(up);
    const double scale = up - lo;

    for (double ratio : linspace(0.0, 1.0, 191)) {
        CHECK_EQ(
            cuberp(start, end, ddstart, ddend, ratio, scale),
            doctest::Approx(func(lo + scale * ratio)).epsilon(kEpsilon)
        );
        CHECK_EQ(
            cuberpd(start, end, ddstart, ddend, ratio, scale),
            doctest::Approx(dfunc(lo + scale * ratio)).epsilon(kEpsilon)
        );
    }
}

} // namespace boyle::math
