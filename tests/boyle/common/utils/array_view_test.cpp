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

TEST_SUITE("ArrayViewTest") {
    TEST_CASE("InitializerListTest") {
        const ArrayView<double> const_array_view{231.8, -34.34, 0.28, 10004.0, 0.005764, 0.0};
        ArrayView<double> array_view{231.8, -34.34, 0.28, 10004.0, 0.005764, 0.0};

        REQUIRE_EQ(const_array_view.size(), 6);

        REQUIRE_EQ(array_view.size(), 6);

        CHECK_EQ(const_array_view[0], 231.8);
        CHECK_EQ(const_array_view[1], -34.34);
        CHECK_EQ(const_array_view[2], 0.28);
        CHECK_EQ(const_array_view[3], 10004.0);
        CHECK_EQ(const_array_view[4], 0.005764);
        CHECK_EQ(const_array_view[5], 0.0);

        CHECK_EQ(array_view[0], 231.8);
        CHECK_EQ(array_view[1], -34.34);
        CHECK_EQ(array_view[2], 0.28);
        CHECK_EQ(array_view[3], 10004.0);
        CHECK_EQ(array_view[4], 0.005764);
        CHECK_EQ(array_view[5], 0.0);

        array_view[3] = -0.000455;

        CHECK_NE(array_view[3], 10004.0);
        CHECK_EQ(array_view[3], -0.000455);
    }

    TEST_CASE("RawArrayTest") {
        // NOLINTNEXTLINE(modernize-avoid-c-arrays)
        double array[6] = {231.8, -34.34, 0.28, 10004.0, 0.005764, 0.0};

        const ArrayView<double> const_array_view{6, array};
        ArrayView<double> array_view{6, array};

        REQUIRE_EQ(const_array_view.size(), 6);

        REQUIRE_EQ(array_view.size(), 6);

        for (std::size_t i = 0; i < 6; ++i) {
            CHECK_EQ(const_array_view[i], array[i]);
            CHECK_EQ(array_view[i], array[i]);
        }

        array_view[3] = -0.000455;

        CHECK_NE(array[3], 10004.0);
        CHECK_EQ(array[3], -0.000455);
        CHECK_NE(const_array_view[3], 10004.0);
        CHECK_EQ(const_array_view[3], -0.000455);
    }

    TEST_CASE_TEMPLATE("StaticArrayTest", T, std::array<double, 6>, Eigen::Vector<double, 6>) {
        T array;
        array[0] = 231.8;
        array[1] = -34.34;
        array[2] = 0.28;
        array[3] = 10004.0;
        array[4] = 0.005764;
        array[5] = 0.0;

        const ArrayView<double> const_array_view{array};
        ArrayView<double> array_view{array};

        REQUIRE_EQ(const_array_view.size(), 6);

        REQUIRE_EQ(array_view.size(), 6);

        for (std::size_t i = 0; i < 6; ++i) {
            CHECK_EQ(const_array_view[i], array[i]);
            CHECK_EQ(array_view[i], array[i]);
        }

        array_view[3] = -0.000455;

        CHECK_NE(array[3], 10004.0);
        CHECK_EQ(array[3], -0.000455);
        CHECK_NE(const_array_view[3], 10004.0);
        CHECK_EQ(const_array_view[3], -0.000455);
    }

    TEST_CASE_TEMPLATE("DynamicArrayTest", T, std::vector<double>, Eigen::Vector<double, Eigen::Dynamic>) {
        T array(6);
        array[0] = 231.8;
        array[1] = -34.34;
        array[2] = 0.28;
        array[3] = 10004.0;
        array[4] = 0.005764;
        array[5] = 0.0;

        const ArrayView<double> const_array_view{array};
        ArrayView<double> array_view{array};

        REQUIRE_EQ(const_array_view.size(), 6);

        REQUIRE_EQ(array_view.size(), 6);

        for (std::size_t i = 0; i < 6; ++i) {
            CHECK_EQ(const_array_view[i], array[i]);
            CHECK_EQ(array_view[i], array[i]);
        }

        array_view[3] = -0.000455;

        CHECK_NE(array[3], 10004.0);
        CHECK_EQ(array[3], -0.000455);
        CHECK_NE(const_array_view[3], 10004.0);
        CHECK_EQ(const_array_view[3], -0.000455);
    }
}

} // namespace boyle::common
