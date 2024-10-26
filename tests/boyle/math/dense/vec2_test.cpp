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

#include "boyle/math/dense/vec2.hpp"

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
    "Constructor", T, Vec2<std::int32_t>, Vec2<std::int64_t>, Vec2s, Vec2d, Vec2c, Vec2z
) {
    const T a(0.0, 0.0);
    CHECK_EQ(a.x, static_cast<typename T::value_type>(0.0));
    CHECK_EQ(a.y, static_cast<typename T::value_type>(0.0));

    const T b(4.0);
    CHECK_EQ(b.x, static_cast<typename T::value_type>(4.0));
    CHECK_EQ(b.y, static_cast<typename T::value_type>(4.0));

    const T c(std::make_pair(1.0, 6.0));
    CHECK_EQ(c.x, static_cast<typename T::value_type>(1.0));
    CHECK_EQ(c.y, static_cast<typename T::value_type>(6.0));

    const T d(c);
    CHECK_EQ(d.x, c.x);
    CHECK_EQ(d.y, c.y);
}

TEST_CASE_TEMPLATE("Basic", T, Vec2s, Vec2d) {
    const T a(1.0, std::numbers::sqrt3);
    CHECK_EQ(a.euclidean(), static_cast<typename T::value_type>(2.0));
    CHECK_EQ(a.angle(), static_cast<typename T::value_type>(M_PI / 3.0));
    CHECK_EQ(std::hypot(a), static_cast<typename T::value_type>(2.0));
    CHECK_EQ(a[0], a.x);
    CHECK_EQ(a[1], a.y);

    const T b = a.rotate(static_cast<typename T::value_type>(M_PI / 6.0));
    CHECK_EQ(b.x, doctest::Approx(0.0).epsilon(1e-6));
    CHECK_EQ(b.y, static_cast<typename T::value_type>(2.0));
    CHECK_EQ(b.angle(), static_cast<typename T::value_type>(M_PI_2));

    const T c = a.normalized();
    CHECK_EQ(c.x, static_cast<typename T::value_type>(0.5));
    CHECK_EQ(c.y, static_cast<typename T::value_type>(std::sqrt(3.0) * 0.5));
}

TEST_CASE_TEMPLATE("IntegralArithmetic", T, Vec2<std::int32_t>, Vec2<std::int64_t>) {
    constexpr T a(94, -37);
    constexpr T b(-18, 27);
    const T result = (a * 10) + (b / 5) - (a * 4);

    CHECK_EQ(result.x, a.x * 6 + b.x / 5);
    CHECK_EQ(result.y, a.y * 6 + b.y / 5);

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y;
    CHECK_EQ(a_dot_b, exact_a_dot_b);
}

TEST_CASE_TEMPLATE("FloatingArithmetic", T, Vec2s, Vec2d) {
    constexpr T a(1.0, 0.0);
    constexpr T b(0.5, std::numbers::sqrt3 * 0.5);
    const T result = (a * 10.5634) + (b / 5.34) - (a * 0.67);

    CHECK_EQ(
        result.x, a.x * static_cast<typename T::value_type>(9.8934) +
                      b.x / static_cast<typename T::value_type>(5.34)
    );
    CHECK_EQ(
        result.y, a.y * static_cast<typename T::value_type>(9.8934) +
                      b.y / static_cast<typename T::value_type>(5.34)
    );

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y;
    CHECK_EQ(a_dot_b, exact_a_dot_b);

    const typename T::value_type a_crossproj_b = a.crossProj(b);
    const typename T::value_type exact_a_crossproj_b = a.x * b.y - a.y * b.x;
    CHECK_EQ(a_crossproj_b, exact_a_crossproj_b);
}

TEST_CASE_TEMPLATE("ComplexArithmetic", T, Vec2c, Vec2z) {
    constexpr T a(typename T::value_type(9.573, 364.26), typename T::value_type(-423.4, 17.44));
    constexpr T b(typename T::value_type(1.545, -5.35), typename T::value_type(49.234, 39.23));
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

    const typename T::value_type a_dot_b = a.dot(b);
    const typename T::value_type exact_a_dot_b = a.x * b.x + a.y * b.y;
    CHECK_EQ(a_dot_b.real(), exact_a_dot_b.real());
    CHECK_EQ(a_dot_b.imag(), exact_a_dot_b.imag());

    const typename T::value_type a_crossproj_b = a.crossProj(b);
    const typename T::value_type exact_a_crossproj_b = a.x * b.y - a.y * b.x;
    CHECK_EQ(a_crossproj_b.real(), exact_a_crossproj_b.real());
    CHECK_EQ(a_crossproj_b.imag(), exact_a_crossproj_b.imag());
}

TEST_CASE("OstreamOperator") {
    constexpr Vec2d a(1274.12, 4454.23);

    std::ostringstream oss;
    oss << a;

    CHECK_EQ(oss.str(), "(x: 1274.12, y: 4454.23)");
}

TEST_CASE("Format") {
    constexpr Vec2d a(1274.12, 4454.23);

    CHECK_EQ(std::format("{}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(std::format("{0}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(std::format("{0:.2}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(std::format("{0:12}", a), "(x:  1274.120000, y:  4454.230000)");
    CHECK_EQ(std::format("{0:12.2}", a), "(x:      1274.12, y:      4454.23)");

    CHECK_EQ(fmt::format("{}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(fmt::format("{0}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(fmt::format("{0:.2}", a), "(x: 1274.12, y: 4454.23)");
    CHECK_EQ(fmt::format("{0:12}", a), "(x:  1274.120000, y:  4454.230000)");
    CHECK_EQ(fmt::format("{0:12.2}", a), "(x:      1274.12, y:      4454.23)");
}

TEST_CASE_TEMPLATE("Serialization", T, Vec2s, Vec2d, Vec2c, Vec2z) {
    constexpr T a(1274.12, 4454.23);

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
