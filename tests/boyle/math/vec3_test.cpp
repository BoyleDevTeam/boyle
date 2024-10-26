/**
 * @file vec3_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-17
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/vec3.hpp"

#include <cmath>
#include <format>
#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "fmt/format.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Constructor") {
    Vec3d a(0.0, 0.0, 0.0);
    CHECK_EQ(a.x, 0.0);
    CHECK_EQ(a.y, 0.0);
    CHECK_EQ(a.z, 0.0);
    CHECK_EQ(sizeof(a.x), 8);
    CHECK_EQ(sizeof(a.y), 8);
    CHECK_EQ(sizeof(a.z), 8);

    Vec3d b(a);
    CHECK_EQ(b.x, a.x);
    CHECK_EQ(b.y, a.y);
    CHECK_EQ(b.z, a.z);

    Vec3f c(a);
    CHECK_EQ(c.x, 0.0F);
    CHECK_EQ(c.y, 0.0F);
    CHECK_EQ(c.z, 0.0F);
    CHECK_EQ(sizeof(c.x), 4);
    CHECK_EQ(sizeof(c.y), 4);
    CHECK_EQ(sizeof(c.z), 4);

    a.x = 3.0;
    a.y = 4.0;
    a.z = 5.0;
    c = b = a;
    CHECK_EQ(c.x, 3.0);
    CHECK_EQ(c.y, 4.0);
    CHECK_EQ(c.z, 5.0);
}

TEST_CASE("Basic") {
    const Vec3d a(1.0, std::sqrt(3.0), std::sqrt(5.0));
    CHECK_EQ(a.euclidean(), 3.0);
    CHECK_EQ(std::hypot(a), 3.0);

    const Vec3d c = a.normalized();
    CHECK_EQ(c.x, 1.0 / 3.0);
    CHECK_EQ(c.y, std::sqrt(3.0) / 3.0);
    CHECK_EQ(c.z, std::sqrt(5.0) / 3.0);
}

TEST_CASE("Arithmetic") {
    const double sqrt_3 = std::sqrt(3.0);
    const Vec3d a(1.0, 0.0, 1.0);
    const Vec3d b(0.5, sqrt_3 * 0.5, -1.0);

    Vec3d c = a + b;
    CHECK_EQ(c.x, 1.5);
    CHECK_EQ(c.y, sqrt_3 * 0.5);
    CHECK_EQ(c.z, 0.0);
    CHECK_EQ(c.euclidean(), sqrt_3);

    Vec3d d = c - b;
    CHECK_EQ(d.x, a.x);
    CHECK_EQ(d.y, a.y);
    CHECK_EQ(d.z, a.z);

    c -= b;
    CHECK_EQ(c.x, a.x);
    CHECK_EQ(c.y, a.y);
    CHECK_EQ(c.z, a.z);

    c += b;
    CHECK_EQ(c.x, 1.5);
    CHECK_EQ(c.y, sqrt_3 * 0.5);
    CHECK_EQ(c.z, 0.0);

    c = a * 0.5;
    CHECK_EQ(c.x, a.x * 0.5);
    CHECK_EQ(c.y, a.y * 0.5);
    CHECK_EQ(c.z, a.z * 0.5);

    c = 0.5 * a;
    CHECK_EQ(c.x, a.x * 0.5);
    CHECK_EQ(c.y, a.y * 0.5);
    CHECK_EQ(c.z, a.z * 0.5);

    c = a;
    c *= 0.5;
    CHECK_EQ(c.x, a.x * 0.5);
    CHECK_EQ(c.y, a.y * 0.5);
    CHECK_EQ(c.z, a.z * 0.5);

    c = a / 2.0;
    CHECK_EQ(c.x, a.x / 2.0);
    CHECK_EQ(c.y, a.y / 2.0);
    CHECK_EQ(c.z, a.z / 2.0);

    c = a;
    c /= 2.0;
    CHECK_EQ(c.x, a.x / 2.0);
    CHECK_EQ(c.y, a.y / 2.0);
    CHECK_EQ(c.z, a.z / 2.0);

    const double a_dot_b = a.dot(b);
    CHECK_EQ(a_dot_b, -0.5);

    const Vec3d b_cross_c = b.cross(c);
    CHECK_EQ(b_cross_c.x, sqrt_3 * 0.25);
    CHECK_EQ(b_cross_c.y, -0.75);
    CHECK_EQ(b_cross_c.z, -sqrt_3 * 0.25);
}

TEST_CASE("OstreamOperator") {
    const Vec3d a(1274.12, 4454.23, 1289.24);

    std::ostringstream oss;
    oss << a;

    CHECK_EQ(oss.str(), "(x: 1274.12, y: 4454.23, z: 1289.24)");
}

TEST_CASE("Format") {
    constexpr Vec3d a(1274.12, 4454.23, -23.5745);

    CHECK_EQ(std::format("{}", a), "(x: 1274.12, y: 4454.23, z: -23.5745)");
    CHECK_EQ(std::format("{0}", a), "(x: 1274.12, y: 4454.23, z: -23.5745)");
    CHECK_EQ(std::format("{0:.2}", a), "(x: 1274.12, y: 4454.23, z: -23.57)");
    CHECK_EQ(std::format("{0:12}", a), "(x:  1274.120000, y:  4454.230000, z:   -23.574500)");
    CHECK_EQ(std::format("{0:12.2}", a), "(x:      1274.12, y:      4454.23, z:       -23.57)");

    CHECK_EQ(fmt::format("{}", a), "(x: 1274.12, y: 4454.23, z: -23.5745)");
    CHECK_EQ(fmt::format("{0}", a), "(x: 1274.12, y: 4454.23, z: -23.5745)");
    CHECK_EQ(fmt::format("{0:.2}", a), "(x: 1274.12, y: 4454.23, z: -23.57)");
    CHECK_EQ(fmt::format("{0:12}", a), "(x:  1274.120000, y:  4454.230000, z:   -23.574500)");
    CHECK_EQ(fmt::format("{0:12.2}", a), "(x:      1274.12, y:      4454.23, z:       -23.57)");
}

TEST_CASE("Serialization") {
    constexpr Vec3d a(1274.12, 4454.23, 1289.24);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << a;

    Vec3d b;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> b;

    CHECK_EQ(a.x, b.x);
    CHECK_EQ(a.y, b.y);
}

} // namespace boyle::math
