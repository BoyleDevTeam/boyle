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

#include "boyle/math/dense/vec3.hpp"

#include <cmath>
#include <format>
#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "fmt/format.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE_TEMPLATE(
    "Constructor", T, Vec3<std::int32_t>, Vec3<std::int64_t>, Vec3s, Vec3d, Vec3s, Vec3z
) {
    const T a(0.0, 0.0, 0.0);
    CHECK_EQ(a.x, static_cast<typename T::value_type>(0.0));
    CHECK_EQ(a.y, static_cast<typename T::value_type>(0.0));
    CHECK_EQ(a.z, static_cast<typename T::value_type>(0.0));

    const T b(4.0);
    CHECK_EQ(b.x, static_cast<typename T::value_type>(4.0));
    CHECK_EQ(b.y, static_cast<typename T::value_type>(4.0));
    CHECK_EQ(b.z, static_cast<typename T::value_type>(4.0));

    const T c(std::make_tuple(1.0, 6.0, 9.0));
    CHECK_EQ(c.x, static_cast<typename T::value_type>(1.0));
    CHECK_EQ(c.y, static_cast<typename T::value_type>(6.0));
    CHECK_EQ(c.z, static_cast<typename T::value_type>(9.0));

    const T d(c);
    CHECK_EQ(d.x, c.x);
    CHECK_EQ(d.y, c.y);
    CHECK_EQ(d.z, c.z);
}

TEST_CASE_TEMPLATE("Basic", T, Vec3s, Vec3d) {
    const T a(1.0, std::sqrt(3.0), std::sqrt(5.0));
    CHECK_EQ(a.euclidean(), doctest::Approx(3.0).epsilon(2E-7));
    CHECK_EQ(std::hypot(a), doctest::Approx(3.0).epsilon(2E-7));
    CHECK_EQ(a[0], a.x);
    CHECK_EQ(a[1], a.y);
    CHECK_EQ(a[2], a.z);

    const T b = a.normalized();
    CHECK_EQ(b.x, doctest::Approx(1.0 / 3.0).epsilon(2E-7));
    CHECK_EQ(b.y, doctest::Approx(std::sqrt(3.0) / 3.0).epsilon(2E-7));
    CHECK_EQ(b.z, doctest::Approx(std::sqrt(5.0) / 3.0).epsilon(2E-7));
}

TEST_CASE_TEMPLATE("IntegralArithmetic", T, Vec3<std::int32_t>, Vec3<std::int64_t>) {
    constexpr T a(94, -37, 42);
    constexpr T b(-18, 27, 75);
    const T result = (a * 10) + (b / 5) - (a * 4);

    CHECK_EQ(result.x, a.x * 6 + b.x / 5);
    CHECK_EQ(result.y, a.y * 6 + b.y / 5);
    CHECK_EQ(result.z, a.z * 6 + b.z / 5);

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y + a.z * b.z;
    CHECK_EQ(a_dot_b, exact_a_dot_b);
}

TEST_CASE_TEMPLATE("FloatingArithmetic", T, Vec3s, Vec3d) {
    constexpr T a(1.0, 0.0, 2.13);
    constexpr T b(0.5, std::numbers::sqrt3 * 0.5, 1.45);
    const T result = (a * 10.5634) + (b / 5.34) - (a * 0.67);

    CHECK_EQ(
        result.x, a.x * static_cast<typename T::value_type>(9.8934) +
                      b.x / static_cast<typename T::value_type>(5.34)
    );
    CHECK_EQ(
        result.y, a.y * static_cast<typename T::value_type>(9.8934) +
                      b.y / static_cast<typename T::value_type>(5.34)
    );
    CHECK_EQ(
        result.z, doctest::Approx(
                      a.z * static_cast<typename T::value_type>(9.8934) +
                      b.z / static_cast<typename T::value_type>(5.34)
                  )
                      .epsilon(2E-7)
    );

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y + a.z * b.z;
    CHECK_EQ(a_dot_b, exact_a_dot_b);

    const T a_cross_b = a.cross(b);
    const typename T::value_type exact_a_cross_b_x = a.y * b.z - a.z * b.y;
    const typename T::value_type exact_a_cross_b_y = a.z * b.x - a.x * b.z;
    const typename T::value_type exact_a_cross_b_z = a.x * b.y - a.y * b.x;
    CHECK_EQ(a_cross_b.x, exact_a_cross_b_x);
    CHECK_EQ(a_cross_b.y, exact_a_cross_b_y);
    CHECK_EQ(a_cross_b.z, exact_a_cross_b_z);

    const typename T::value_type a_crossproj_b = a.crossProj(b);
    const typename T::value_type exact_a_crossproj_b = a_cross_b.euclidean();
    CHECK_EQ(a_crossproj_b, exact_a_crossproj_b);
}

TEST_CASE_TEMPLATE("ComplexArithmetic", T, Vec3c, Vec3z) {
    constexpr T a(
        typename T::value_type(9.573, 364.26), typename T::value_type(-423.4, 17.44),
        typename T::value_type(23.635, 0.464)
    );
    constexpr T b(
        typename T::value_type(1.545, -5.35), typename T::value_type(49.234, 39.23),
        typename T::value_type(-0.564, -4.65)
    );
    const T result = (a * 10.5634) + (b / 5.34) - (a * 0.67);

    CHECK_EQ(
        result.x.real(), a.x.real() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.x.real() / static_cast<typename T::value_type::value_type>(5.34)
    );
    CHECK_EQ(
        result.x.imag(), a.x.imag() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.x.imag() / static_cast<typename T::value_type::value_type>(5.34)
    );
    CHECK_EQ(
        result.y.real(), doctest::Approx(
                             a.y.real() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.y.real() / static_cast<typename T::value_type::value_type>(5.34)
                         )
                             .epsilon(2E-7)
    );
    CHECK_EQ(
        result.y.imag(), a.y.imag() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.y.imag() / static_cast<typename T::value_type::value_type>(5.34)
    );
    CHECK_EQ(
        result.z.real(), doctest::Approx(
                             a.z.real() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.z.real() / static_cast<typename T::value_type::value_type>(5.34)
                         )
                             .epsilon(2E-7)
    );
    CHECK_EQ(
        result.z.imag(), a.z.imag() * static_cast<typename T::value_type::value_type>(9.8934) +
                             b.z.imag() / static_cast<typename T::value_type::value_type>(5.34)
    );

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y + a.z * b.z;
    CHECK_EQ(a_dot_b.real(), exact_a_dot_b.real());
    CHECK_EQ(a_dot_b.imag(), exact_a_dot_b.imag());

    const T a_cross_b = a.cross(b);
    const typename T::value_type exact_a_cross_b_x = a.y * b.z - a.z * b.y;
    const typename T::value_type exact_a_cross_b_y = a.z * b.x - a.x * b.z;
    const typename T::value_type exact_a_cross_b_z = a.x * b.y - a.y * b.x;
    CHECK_EQ(a_cross_b.x.real(), exact_a_cross_b_x.real());
    CHECK_EQ(a_cross_b.x.imag(), exact_a_cross_b_x.imag());
    CHECK_EQ(a_cross_b.y.real(), exact_a_cross_b_y.real());
    CHECK_EQ(a_cross_b.y.imag(), exact_a_cross_b_y.imag());
    CHECK_EQ(a_cross_b.z.real(), exact_a_cross_b_z.real());
    CHECK_EQ(a_cross_b.z.imag(), exact_a_cross_b_z.imag());

    const typename T::value_type a_crossproj_b = a.crossProj(b);
    const typename T::value_type exact_a_crossproj_b = a_cross_b.euclidean();
    CHECK_EQ(a_crossproj_b.real(), exact_a_crossproj_b.real());
    CHECK_EQ(a_crossproj_b.imag(), exact_a_crossproj_b.imag());
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

TEST_CASE_TEMPLATE("Serialization", T, Vec3s, Vec3d, Vec3c, Vec3z) {
    constexpr T a(1274.12, 4454.23, 1289.24);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << a;

    T b;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> b;

    CHECK_EQ(a.x, b.x);
    CHECK_EQ(a.y, b.y);
}

} // namespace boyle::math
