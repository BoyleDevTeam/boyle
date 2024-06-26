/**
 * @file vec2_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-17
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/vec2.hpp"

#include <cmath>
#include <sstream>
#include <string>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Constructor") {
    Vec2d a(0.0, 0.0);
    CHECK_EQ(a.x(), 0.0);
    CHECK_EQ(a.y(), 0.0);
    CHECK_EQ(sizeof(a.x()), 8);
    CHECK_EQ(sizeof(a.y()), 8);

    Vec2d b(a);
    CHECK_EQ(b.x(), a.x());
    CHECK_EQ(b.y(), a.y());

    Vec2f c(a);
    CHECK_EQ(c.x(), 0.0f);
    CHECK_EQ(c.y(), 0.0f);
    CHECK_EQ(sizeof(c.x()), 4);
    CHECK_EQ(sizeof(c.y()), 4);

    a.x() = 3.0;
    a.y() = 4.0;
    c = b = a;
    CHECK_EQ(c.x(), 3.0);
    CHECK_EQ(c.y(), 4.0);
}

TEST_CASE("Basic") {
    const Vec2d a(1.0, std::sqrt(3.0));
    CHECK_EQ(a.norm(), 2.0);
    CHECK_EQ(a.angle(), M_PI / 3.0);
    CHECK_EQ(std::hypot(a), 2.0);

    const Vec2d b = a.rotate(M_PI / 6.0);
    CHECK_EQ(b.x(), doctest::Approx(0.0).epsilon(1e-6));
    CHECK_EQ(b.y(), 2.0);
    CHECK_EQ(b.angle(), M_PI_2);

    const Vec2d c = normalize(a);
    CHECK_EQ(c.x(), 0.5);
    CHECK_EQ(c.y(), std::sqrt(3.0) * 0.5);
}

TEST_CASE("Arithmetic") {
    const double sqrt_3 = std::sqrt(3.0);
    const Vec2d a(1.0, 0.0);
    const Vec2d b(0.5, sqrt_3 * 0.5);

    Vec2d c = a + b;
    CHECK_EQ(c.x(), 1.5);
    CHECK_EQ(c.y(), sqrt_3 * 0.5);
    CHECK_EQ(c.norm(), sqrt_3);
    CHECK_EQ(c.angle(), M_PI / 6.0);

    Vec2d d = c - b;
    CHECK_EQ(d.x(), a.x());
    CHECK_EQ(d.y(), a.y());

    c -= b;
    CHECK_EQ(c.x(), a.x());
    CHECK_EQ(c.y(), a.y());

    c += b;
    CHECK_EQ(c.x(), 1.5);
    CHECK_EQ(c.y(), sqrt_3 * 0.5);

    c = a * 0.5;
    CHECK_EQ(c.x(), a.x() * 0.5);
    CHECK_EQ(c.y(), a.y() * 0.5);

    c = 0.5 * a;
    CHECK_EQ(c.x(), a.x() * 0.5);
    CHECK_EQ(c.y(), a.y() * 0.5);

    c = a;
    c *= 0.5;
    CHECK_EQ(c.x(), a.x() * 0.5);
    CHECK_EQ(c.y(), a.y() * 0.5);

    c = a / 2.0;
    CHECK_EQ(c.x(), a.x() / 2.0);
    CHECK_EQ(c.y(), a.y() / 2.0);

    c = a;
    c /= 2.0;
    CHECK_EQ(c.x(), a.x() / 2.0);
    CHECK_EQ(c.y(), a.y() / 2.0);

    const double a_dot_b = a.dot(b);
    CHECK_EQ(a_dot_b, 0.5);

    const double b_cross_c = b.cross(c);
    CHECK_EQ(b_cross_c, -sqrt_3 * 0.25);
}

TEST_CASE("OstreamOperator") {
    const Vec2d a(1274.12, 4454.23);

    std::ostringstream oss;
    oss << a;

    const std::string a_text = "(1274.12, 4454.23)";

    CHECK_EQ(oss.str(), a_text);
}

TEST_CASE("Serialization") {
    const Vec2d a(1274.12, 4454.23);
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << a;

    Vec2d b;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> b;

    CHECK_EQ(a.x(), b.x());
    CHECK_EQ(a.y(), b.y());
}

} // namespace boyle::math
