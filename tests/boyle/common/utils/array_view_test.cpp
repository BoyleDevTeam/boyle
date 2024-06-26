/**
 * @file array_view_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/common/utils/array_view.hpp"

#include <array>
#include <vector>

#include "Eigen/Eigen"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "boyle/common/utils/macros.hpp"

namespace boyle::common {

class ArrayViewTestFixture {
  public:
    ArrayViewTestFixture() noexcept : std_vector_arr(6), eigen_dynamic_vector_arr(6) {
        raw_arr[0] = 231.8;
        raw_arr[1] = -34.34;
        raw_arr[2] = 0.28;
        raw_arr[3] = 10004.0;
        raw_arr[4] = 0.005764;
        raw_arr[5] = 0.0;
        std_array_arr[0] = 231.8;
        std_array_arr[1] = -34.34;
        std_array_arr[2] = 0.28;
        std_array_arr[3] = 10004.0;
        std_array_arr[4] = 0.005764;
        std_array_arr[5] = 0.0;
        std_vector_arr[0] = 231.8;
        std_vector_arr[1] = -34.34;
        std_vector_arr[2] = 0.28;
        std_vector_arr[3] = 10004.0;
        std_vector_arr[4] = 0.005764;
        std_vector_arr[5] = 0.0;
        eigen_static_vector_arr(0) = 231.8;
        eigen_static_vector_arr(1) = -34.34;
        eigen_static_vector_arr(2) = 0.28;
        eigen_static_vector_arr(3) = 10004.0;
        eigen_static_vector_arr(4) = 0.005764;
        eigen_static_vector_arr(5) = 0.0;
        eigen_dynamic_vector_arr(0) = 231.8;
        eigen_dynamic_vector_arr(1) = -34.34;
        eigen_dynamic_vector_arr(2) = 0.28;
        eigen_dynamic_vector_arr(3) = 10004.0;
        eigen_dynamic_vector_arr(4) = 0.005764;
        eigen_dynamic_vector_arr(5) = 0.0;
    }
    DISABLE_COPY_AND_MOVE(ArrayViewTestFixture);
    ~ArrayViewTestFixture() noexcept = default;

  protected:
    double raw_arr[6]; // NOLINT(modernize-avoid-c-arrays)
    std::array<double, 6> std_array_arr;
    std::vector<double> std_vector_arr;
    Eigen::Vector<double, 6> eigen_static_vector_arr;
    Eigen::VectorXd eigen_dynamic_vector_arr;
};

TEST_SUITE("DoubleArrayTraverseRead") {
    TEST_CASE_FIXTURE(ArrayViewTestFixture, "Constructor") {
        SUBCASE("std_initialize_list") {
            const ArrayView<double> array_view{231.8, -34.34, 0.28, 10004.0, 0.005764, 0.0};
            REQUIRE_EQ(array_view.size(), 6);
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], raw_arr[i]);
            }
        }
        SUBCASE("raw_array_1") {
            const ArrayView<double> array_view(raw_arr, 6);
            REQUIRE_EQ(array_view.data(), raw_arr);
            REQUIRE_EQ(array_view.size(), 6);
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], raw_arr[i]);
            }
        }
        SUBCASE("raw_array_2") {
            const ArrayView<double> array_view(raw_arr, raw_arr + 6);
            REQUIRE_EQ(array_view.data(), raw_arr);
            REQUIRE_EQ(array_view.size(), 6);
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], raw_arr[i]);
            }
        }
        SUBCASE("std_array") {
            const ArrayView<double> array_view(std_array_arr);
            REQUIRE_EQ(array_view.data(), std_array_arr.data());
            REQUIRE_EQ(array_view.size(), std_array_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], std_array_arr[i]);
            }
        }
        SUBCASE("std_vector") {
            const ArrayView<double> array_view(std_vector_arr);
            REQUIRE_EQ(array_view.data(), std_vector_arr.data());
            REQUIRE_EQ(array_view.size(), std_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], std_vector_arr[i]);
            }
        }
        SUBCASE("eigen_static_vector") {
            const ArrayView<double> array_view(eigen_static_vector_arr);
            REQUIRE_EQ(array_view.data(), eigen_static_vector_arr.data());
            REQUIRE_EQ(array_view.size(), eigen_static_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], eigen_static_vector_arr(i));
            }
        }
        SUBCASE("eigen_dynamic_vector") {
            const ArrayView<double> array_view(eigen_dynamic_vector_arr);
            REQUIRE_EQ(array_view.data(), eigen_dynamic_vector_arr.data());
            REQUIRE_EQ(array_view.size(), eigen_dynamic_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], eigen_dynamic_vector_arr(i));
            }
        }
    }
}

TEST_SUITE("DoubleArrayTraverseWrite") {
    TEST_CASE_FIXTURE(ArrayViewTestFixture, "Constructor") {
        SUBCASE("raw_array_1") {
            ArrayView<double> array_view(raw_arr, 6);
            REQUIRE_EQ(array_view.data(), raw_arr);
            REQUIRE_EQ(array_view.size(), 6);
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], raw_arr[i]);
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], raw_arr[3]);
        }
        SUBCASE("raw_array_2") {
            ArrayView<double> array_view(raw_arr, raw_arr + 6);
            REQUIRE_EQ(array_view.data(), raw_arr);
            REQUIRE_EQ(array_view.size(), 6);
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], raw_arr[i]);
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], raw_arr[3]);
        }
        SUBCASE("std_array") {
            ArrayView<double> array_view(std_array_arr);
            REQUIRE_EQ(array_view.data(), std_array_arr.data());
            REQUIRE_EQ(array_view.size(), std_array_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], std_array_arr[i]);
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], std_array_arr[3]);
        }
        SUBCASE("std_vector") {
            ArrayView<double> array_view(std_vector_arr);
            REQUIRE_EQ(array_view.data(), std_vector_arr.data());
            REQUIRE_EQ(array_view.size(), std_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], std_vector_arr[i]);
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], std_vector_arr[3]);
        }
        SUBCASE("eigen_static_vector") {
            ArrayView<double> array_view(eigen_static_vector_arr);
            REQUIRE_EQ(array_view.data(), eigen_static_vector_arr.data());
            REQUIRE_EQ(array_view.size(), eigen_static_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], eigen_static_vector_arr(i));
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], eigen_static_vector_arr[3]);
        }
        SUBCASE("eigen_dynamic_vector") {
            ArrayView<double> array_view(eigen_dynamic_vector_arr);
            REQUIRE_EQ(array_view.data(), eigen_dynamic_vector_arr.data());
            REQUIRE_EQ(array_view.size(), eigen_dynamic_vector_arr.size());
            for (std::size_t i = 0; i < array_view.size(); ++i) {
                CHECK_EQ(array_view[i], eigen_dynamic_vector_arr(i));
            }
            array_view[3] = -0.000455;
            CHECK_EQ(array_view[3], eigen_dynamic_vector_arr[3]);
        }
    }
}

} // namespace boyle::common
