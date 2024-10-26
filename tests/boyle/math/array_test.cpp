/**
 * @file array_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-14
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/math/array.hpp"

#include <array>
#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

class ArrayTestFixture {
  public:
    ArrayTestFixture(const ArrayTestFixture& other) noexcept = delete;
    auto operator=(const ArrayTestFixture& other) noexcept -> ArrayTestFixture& = delete;
    ArrayTestFixture(ArrayTestFixture&& other) noexcept = delete;
    auto operator=(ArrayTestFixture&& other) noexcept -> ArrayTestFixture& = delete;
    ~ArrayTestFixture() = default;

    ArrayTestFixture() noexcept
        : double_math_arrays{Array<double>{99.5124, -72.4566, -13.7342, 12.6311, -79.6078, -51.7363}, Array<double>(6, 0.0)},
          float_math_arrays{double_math_arrays[0], Array<float>(6)} {
        float_math_arrays[1][0] = -39.1180F;
        float_math_arrays[1][1] = 38.5798F;
        float_math_arrays[1][2] = 91.1946F;
        float_math_arrays[1][3] = -41.6234F;
        float_math_arrays[1][4] = 87.4647F;
        float_math_arrays[1][5] = 71.3192F;
        double_math_arrays[1] = float_math_arrays[1];
    }

  protected:
    std::array<Array<double>, 2> double_math_arrays;
    std::array<Array<float>, 2> float_math_arrays;
};

TEST_CASE_FIXTURE(ArrayTestFixture, "PlusTest") {
    const Array<double> double_plus_double{double_math_arrays[0] + double_math_arrays[1]};
    CHECK_EQ(double_plus_double[0], doctest::Approx(99.5124 - 39.1180).epsilon(1E-5));
    CHECK_EQ(double_plus_double[1], doctest::Approx(-72.4566 + 38.5798).epsilon(1E-5));
    CHECK_EQ(double_plus_double[2], doctest::Approx(-13.7342 + 91.1946).epsilon(1E-5));
    CHECK_EQ(double_plus_double[3], doctest::Approx(12.6311 - 41.6234).epsilon(1E-5));
    CHECK_EQ(double_plus_double[4], doctest::Approx(-79.6078 + 87.4647).epsilon(1E-5));
    CHECK_EQ(double_plus_double[5], doctest::Approx(-51.7363 + 71.3192).epsilon(1E-5));

    const Array<float> float_plus_float{float_math_arrays[0] + float_math_arrays[1]};
    CHECK_EQ(float_plus_float[0], doctest::Approx(99.5124 - 39.1180).epsilon(1E-5));
    CHECK_EQ(float_plus_float[1], doctest::Approx(-72.4566 + 38.5798).epsilon(1E-5));
    CHECK_EQ(float_plus_float[2], doctest::Approx(-13.7342 + 91.1946).epsilon(1E-5));
    CHECK_EQ(float_plus_float[3], doctest::Approx(12.6311 - 41.6234).epsilon(1E-5));
    CHECK_EQ(float_plus_float[4], doctest::Approx(-79.6078 + 87.4647).epsilon(1E-5));
    CHECK_EQ(float_plus_float[5], doctest::Approx(-51.7363 + 71.3192).epsilon(1E-5));

    const Array<double> double_plus_float{double_math_arrays[0] + float_math_arrays[1]};
    CHECK_EQ(double_plus_float[0], doctest::Approx(99.5124 - 39.1180).epsilon(1E-5));
    CHECK_EQ(double_plus_float[1], doctest::Approx(-72.4566 + 38.5798).epsilon(1E-5));
    CHECK_EQ(double_plus_float[2], doctest::Approx(-13.7342 + 91.1946).epsilon(1E-5));
    CHECK_EQ(double_plus_float[3], doctest::Approx(12.6311 - 41.6234).epsilon(1E-5));
    CHECK_EQ(double_plus_float[4], doctest::Approx(-79.6078 + 87.4647).epsilon(1E-5));
    CHECK_EQ(double_plus_float[5], doctest::Approx(-51.7363 + 71.3192).epsilon(1E-5));

    const Array<double> float_plus_double{float_math_arrays[0] + double_math_arrays[1]};
    CHECK_EQ(float_plus_double[0], doctest::Approx(99.5124 - 39.1180).epsilon(1E-5));
    CHECK_EQ(float_plus_double[1], doctest::Approx(-72.4566 + 38.5798).epsilon(1E-5));
    CHECK_EQ(float_plus_double[2], doctest::Approx(-13.7342 + 91.1946).epsilon(1E-5));
    CHECK_EQ(float_plus_double[3], doctest::Approx(12.6311 - 41.6234).epsilon(1E-5));
    CHECK_EQ(float_plus_double[4], doctest::Approx(-79.6078 + 87.4647).epsilon(1E-5));
    CHECK_EQ(float_plus_double[5], doctest::Approx(-51.7363 + 71.3192).epsilon(1E-5));
}

TEST_CASE_FIXTURE(ArrayTestFixture, "MinusTest") {
    const Array<double> double_minus_double{double_math_arrays[0] - double_math_arrays[1]};
    CHECK_EQ(double_minus_double[0], doctest::Approx(99.5124 + 39.1180).epsilon(1E-5));
    CHECK_EQ(double_minus_double[1], doctest::Approx(-72.4566 - 38.5798).epsilon(1E-5));
    CHECK_EQ(double_minus_double[2], doctest::Approx(-13.7342 - 91.1946).epsilon(1E-5));
    CHECK_EQ(double_minus_double[3], doctest::Approx(12.6311 + 41.6234).epsilon(1E-5));
    CHECK_EQ(double_minus_double[4], doctest::Approx(-79.6078 - 87.4647).epsilon(1E-5));
    CHECK_EQ(double_minus_double[5], doctest::Approx(-51.7363 - 71.3192).epsilon(1E-5));

    const Array<float> float_minus_float{float_math_arrays[0] - float_math_arrays[1]};
    CHECK_EQ(float_minus_float[0], doctest::Approx(99.5124 + 39.1180).epsilon(1E-5));
    CHECK_EQ(float_minus_float[1], doctest::Approx(-72.4566 - 38.5798).epsilon(1E-5));
    CHECK_EQ(float_minus_float[2], doctest::Approx(-13.7342 - 91.1946).epsilon(1E-5));
    CHECK_EQ(float_minus_float[3], doctest::Approx(12.6311 + 41.6234).epsilon(1E-5));
    CHECK_EQ(float_minus_float[4], doctest::Approx(-79.6078 - 87.4647).epsilon(1E-5));
    CHECK_EQ(float_minus_float[5], doctest::Approx(-51.7363 - 71.3192).epsilon(1E-5));

    const Array<double> double_minus_float{double_math_arrays[0] - float_math_arrays[1]};
    CHECK_EQ(double_minus_float[0], doctest::Approx(99.5124 + 39.1180).epsilon(1E-5));
    CHECK_EQ(double_minus_float[1], doctest::Approx(-72.4566 - 38.5798).epsilon(1E-5));
    CHECK_EQ(double_minus_float[2], doctest::Approx(-13.7342 - 91.1946).epsilon(1E-5));
    CHECK_EQ(double_minus_float[3], doctest::Approx(12.6311 + 41.6234).epsilon(1E-5));
    CHECK_EQ(double_minus_float[4], doctest::Approx(-79.6078 - 87.4647).epsilon(1E-5));
    CHECK_EQ(double_minus_float[5], doctest::Approx(-51.7363 - 71.3192).epsilon(1E-5));

    const Array<double> float_minus_double{float_math_arrays[0] - double_math_arrays[1]};
    CHECK_EQ(float_minus_double[0], doctest::Approx(99.5124 + 39.1180).epsilon(1E-5));
    CHECK_EQ(float_minus_double[1], doctest::Approx(-72.4566 - 38.5798).epsilon(1E-5));
    CHECK_EQ(float_minus_double[2], doctest::Approx(-13.7342 - 91.1946).epsilon(1E-5));
    CHECK_EQ(float_minus_double[3], doctest::Approx(12.6311 + 41.6234).epsilon(1E-5));
    CHECK_EQ(float_minus_double[4], doctest::Approx(-79.6078 - 87.4647).epsilon(1E-5));
    CHECK_EQ(float_minus_double[5], doctest::Approx(-51.7363 - 71.3192).epsilon(1E-5));
}

TEST_CASE_FIXTURE(ArrayTestFixture, "MultiplyTest") {
    const Array<double> double_multiply_double{double_math_arrays[0] * 25.5143};
    CHECK_EQ(double_multiply_double[0], doctest::Approx(99.5124 * 25.5143).epsilon(1E-12));
    CHECK_EQ(double_multiply_double[1], doctest::Approx(-72.4566 * 25.5143).epsilon(1E-12));
    CHECK_EQ(double_multiply_double[2], doctest::Approx(-13.7342 * 25.5143).epsilon(1E-12));
    CHECK_EQ(double_multiply_double[3], doctest::Approx(12.6311 * 25.5143).epsilon(1E-12));
    CHECK_EQ(double_multiply_double[4], doctest::Approx(-79.6078 * 25.5143).epsilon(1E-12));
    CHECK_EQ(double_multiply_double[5], doctest::Approx(-51.7363 * 25.5143).epsilon(1E-12));

    const Array<float> float_multiply_float{float_math_arrays[0] * 25.5143F};
    CHECK_EQ(float_multiply_float[0], doctest::Approx(99.5124F * 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_multiply_float[1], doctest::Approx(-72.4566F * 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_multiply_float[2], doctest::Approx(-13.7342F * 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_multiply_float[3], doctest::Approx(12.6311F * 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_multiply_float[4], doctest::Approx(-79.6078F * 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_multiply_float[5], doctest::Approx(-51.7363F * 25.5143F).epsilon(1E-12));

    const Array<double> double_multiply_float{25.5148F * double_math_arrays[0]};
    CHECK_EQ(double_multiply_float[0], doctest::Approx(99.5124 * 25.5148F).epsilon(1E-12));
    CHECK_EQ(double_multiply_float[1], doctest::Approx(-72.4566 * 25.5148F).epsilon(1E-12));
    CHECK_EQ(double_multiply_float[2], doctest::Approx(-13.7342 * 25.5148F).epsilon(1E-12));
    CHECK_EQ(double_multiply_float[3], doctest::Approx(12.6311 * 25.5148F).epsilon(1E-12));
    CHECK_EQ(double_multiply_float[4], doctest::Approx(-79.6078 * 25.5148F).epsilon(1E-12));
    CHECK_EQ(double_multiply_float[5], doctest::Approx(-51.7363 * 25.5148F).epsilon(1E-12));

    const Array<float> float_multiply_double{25.5148 * float_math_arrays[0]};
    CHECK_EQ(float_multiply_double[0], doctest::Approx(99.5124F * 25.5148).epsilon(1E-6));
    CHECK_EQ(float_multiply_double[1], doctest::Approx(-72.4566F * 25.5148).epsilon(1E-6));
    CHECK_EQ(float_multiply_double[2], doctest::Approx(-13.7342F * 25.5148).epsilon(1E-6));
    CHECK_EQ(float_multiply_double[3], doctest::Approx(12.6311F * 25.5148).epsilon(1E-6));
    CHECK_EQ(float_multiply_double[4], doctest::Approx(-79.6078F * 25.5148).epsilon(1E-6));
    CHECK_EQ(float_multiply_double[5], doctest::Approx(-51.7363F * 25.5148).epsilon(1E-6));
}

TEST_CASE_FIXTURE(ArrayTestFixture, "DivideTest") {
    const Array<double> double_divide_double{double_math_arrays[0] / 25.5143};
    CHECK_EQ(double_divide_double[0], doctest::Approx(99.5124 / 25.5143).epsilon(1E-12));
    CHECK_EQ(double_divide_double[1], doctest::Approx(-72.4566 / 25.5143).epsilon(1E-12));
    CHECK_EQ(double_divide_double[2], doctest::Approx(-13.7342 / 25.5143).epsilon(1E-12));
    CHECK_EQ(double_divide_double[3], doctest::Approx(12.6311 / 25.5143).epsilon(1E-12));
    CHECK_EQ(double_divide_double[4], doctest::Approx(-79.6078 / 25.5143).epsilon(1E-12));
    CHECK_EQ(double_divide_double[5], doctest::Approx(-51.7363 / 25.5143).epsilon(1E-12));

    const Array<float> float_divide_float{float_math_arrays[0] / 25.5143F};
    CHECK_EQ(float_divide_float[0], doctest::Approx(99.5124F / 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_divide_float[1], doctest::Approx(-72.4566F / 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_divide_float[2], doctest::Approx(-13.7342F / 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_divide_float[3], doctest::Approx(12.6311F / 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_divide_float[4], doctest::Approx(-79.6078F / 25.5143F).epsilon(1E-12));
    CHECK_EQ(float_divide_float[5], doctest::Approx(-51.7363F / 25.5143F).epsilon(1E-6));

    const Array<double> double_divide_float{double_math_arrays[0] / 25.5143F};
    CHECK_EQ(double_divide_float[0], doctest::Approx(99.5124 / 25.5143F).epsilon(1E-12));
    CHECK_EQ(double_divide_float[1], doctest::Approx(-72.4566 / 25.5143F).epsilon(1E-12));
    CHECK_EQ(double_divide_float[2], doctest::Approx(-13.7342 / 25.5143F).epsilon(1E-12));
    CHECK_EQ(double_divide_float[3], doctest::Approx(12.6311 / 25.5143F).epsilon(1E-12));
    CHECK_EQ(double_divide_float[4], doctest::Approx(-79.6078 / 25.5143F).epsilon(1E-12));
    CHECK_EQ(double_divide_float[5], doctest::Approx(-51.7363 / 25.5143F).epsilon(1E-12));

    const Array<float> float_divide_double{float_math_arrays[0] / 25.5143};
    CHECK_EQ(float_divide_double[0], doctest::Approx(99.5124F / 25.5143).epsilon(1E-6));
    CHECK_EQ(float_divide_double[1], doctest::Approx(-72.4566F / 25.5143).epsilon(1E-6));
    CHECK_EQ(float_divide_double[2], doctest::Approx(-13.7342F / 25.5143).epsilon(1E-6));
    CHECK_EQ(float_divide_double[3], doctest::Approx(12.6311F / 25.5143).epsilon(1E-6));
    CHECK_EQ(float_divide_double[4], doctest::Approx(-79.6078F / 25.5143).epsilon(1E-6));
    CHECK_EQ(float_divide_double[5], doctest::Approx(-51.7363F / 25.5143).epsilon(1E-6));
}

TEST_CASE_FIXTURE(ArrayTestFixture, "DotTest") {
    constexpr double exact_result = 99.5124 * -39.1180 + -72.4566 * 38.5798 + -13.7342 * 91.1946 +
                                    12.6311 * -41.6234 + -79.6078 * 87.4647 + -51.7363 * 71.3192;

    const double double_dot_double{double_math_arrays[0].dot(double_math_arrays[1])};
    CHECK_EQ(double_dot_double, doctest::Approx(exact_result).epsilon(1E-8));

    const float float_dot_float{float_math_arrays[0].dot(float_math_arrays[1])};
    CHECK_EQ(float_dot_float, doctest::Approx(exact_result).epsilon(1E-6));

    const double double_dot_float{double_math_arrays[0].dot(float_math_arrays[1])};
    CHECK_EQ(double_dot_float, doctest::Approx(exact_result).epsilon(1E-8));

    const double float_dot_double{float_math_arrays[0].dot(double_math_arrays[1])};
    CHECK_EQ(float_dot_double, doctest::Approx(exact_result).epsilon(1E-6));
}

TEST_CASE_FIXTURE(ArrayTestFixture, "Serialization") {

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << double_math_arrays[0];
    oa << float_math_arrays[0];

    Array<double> other_double_math_array;
    Array<float> other_float_math_array;

    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_double_math_array;
    ia >> other_float_math_array;

    REQUIRE_EQ(double_math_arrays[0].size(), other_double_math_array.size());

    for (std::size_t i{0}; i < double_math_arrays[0].size(); ++i) {
        CHECK_EQ(double_math_arrays[0][i], other_double_math_array[i]);
    }

    REQUIRE_EQ(float_math_arrays[0].size(), other_float_math_array.size());

    for (std::size_t i{0}; i < float_math_arrays[0].size(); ++i) {
        CHECK_EQ(float_math_arrays[0][i], other_float_math_array[i]);
    }
}

} // namespace boyle::math
