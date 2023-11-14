/**
 * @file csc_matrix_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-11-04
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "math/sparse_matrix/csc_matrix.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "math/sparse_matrix/coo_matrix.hpp"
#include "math/sparse_matrix/lil_matrix.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace tiny_pnc {
namespace math {

TEST_CASE("Basic") {
    CscMatrix csc_matrix(5, 5);

    CHECK_EQ(csc_matrix.nrows(), 5);
    CHECK_EQ(csc_matrix.ncols(), 5);
    CHECK_EQ(csc_matrix.nnzs(), 0);
    CHECK_EQ(csc_matrix.outerIndices().size(), 6);

    csc_matrix.updateCoeff(4, 4, 8.0);
    csc_matrix.updateCoeff(1, 4, 17.0);
    csc_matrix.updateCoeff(2, 3, 1.0);
    csc_matrix.updateCoeff(4, 2, 14.0);
    csc_matrix.updateCoeff(2, 1, 5.0);
    csc_matrix.updateCoeff(0, 1, 3.0);
    csc_matrix.updateCoeff(2, 0, 7.0);
    csc_matrix.updateCoeff(1, 0, 22.0);
    csc_matrix.updateCoeff(5, 5, 1.0);

    CHECK_EQ(csc_matrix.coeff(1, 0), 22.0);
    CHECK_EQ(csc_matrix.coeff(2, 0), 7.0);
    CHECK_EQ(csc_matrix.coeff(0, 1), 3.0);
    CHECK_EQ(csc_matrix.coeff(2, 1), 5.0);
    CHECK_EQ(csc_matrix.coeff(4, 2), 14.0);
    CHECK_EQ(csc_matrix.coeff(2, 3), 1.0);
    CHECK_EQ(csc_matrix.coeff(1, 4), 17.0);
    CHECK_EQ(csc_matrix.coeff(4, 4), 8.0);
    CHECK_EQ(csc_matrix.coeff(5, 5), 0.0);
    CHECK_EQ(csc_matrix.coeff(3, 3), 0.0);

    const std::size_t values_size{csc_matrix.values().size()};
    const std::size_t inner_indices_size{csc_matrix.innerIndices().size()};
    const std::size_t outer_indices_size{csc_matrix.outerIndices().size()};
    const std::vector<double>& values{csc_matrix.values()};
    const std::vector<int>& inner_indices{csc_matrix.innerIndices()};
    const std::vector<int>& outer_indices{csc_matrix.outerIndices()};

    const std::vector<double> exact_values{22.0, 7.0, 3.0, 5.0, 14.0, 1.0, 17.0, 8.0};
    const std::vector<int> exact_inner_indices{1, 2, 0, 2, 4, 2, 1, 4};
    const std::vector<int> exact_outer_indices{0, 2, 4, 5, 6, 8};

    CHECK_EQ(values_size, exact_values.size());
    CHECK_EQ(inner_indices_size, exact_inner_indices.size());
    CHECK_EQ(outer_indices_size, exact_outer_indices.size());

    for (std::size_t i{0}; i < values_size; ++i) {
        CHECK_EQ(values[i], exact_values[i]);
    }
    for (std::size_t i{0}; i < inner_indices_size; ++i) {
        CHECK_EQ(inner_indices[i], exact_inner_indices[i]);
    }
    for (std::size_t i{0}; i < outer_indices_size; ++i) {
        CHECK_EQ(outer_indices[i], exact_outer_indices[i]);
    }
}

TEST_CASE("Conversion") {
    CooMatrix coo_matrix(5, 5);

    coo_matrix.updateCoeff(4, 4, 8.0);
    coo_matrix.updateCoeff(1, 4, 17.0);
    coo_matrix.updateCoeff(2, 3, 1.0);
    coo_matrix.updateCoeff(4, 2, 14.0);
    coo_matrix.updateCoeff(2, 1, 5.0);
    coo_matrix.updateCoeff(0, 1, 3.0);
    coo_matrix.updateCoeff(2, 0, 7.0);
    coo_matrix.updateCoeff(1, 0, 22.0);
    coo_matrix.updateCoeff(5, 5, 1.0);

    CscMatrix csc_matrix{coo_matrix};

    CHECK_EQ(csc_matrix.nrows(), 5);
    CHECK_EQ(csc_matrix.ncols(), 5);
    CHECK_EQ(csc_matrix.nnzs(), 8);
    CHECK_EQ(csc_matrix.outerIndices().size(), 6);

    const std::size_t values_size{csc_matrix.values().size()};
    const std::size_t inner_indices_size{csc_matrix.innerIndices().size()};
    const std::size_t outer_indices_size{csc_matrix.outerIndices().size()};
    const std::vector<double>& values{csc_matrix.values()};
    const std::vector<int>& inner_indices{csc_matrix.innerIndices()};
    const std::vector<int>& outer_indices{csc_matrix.outerIndices()};

    const std::vector<double> exact_values{22.0, 7.0, 3.0, 5.0, 14.0, 1.0, 17.0, 8.0};
    const std::vector<int> exact_inner_indices{1, 2, 0, 2, 4, 2, 1, 4};
    const std::vector<int> exact_outer_indices{0, 2, 4, 5, 6, 8};

    CHECK_EQ(values_size, exact_values.size());
    CHECK_EQ(inner_indices_size, exact_inner_indices.size());
    CHECK_EQ(outer_indices_size, exact_outer_indices.size());

    for (std::size_t i{0}; i < values_size; ++i) {
        CHECK_EQ(values[i], exact_values[i]);
    }
    for (std::size_t i{0}; i < inner_indices_size; ++i) {
        CHECK_EQ(inner_indices[i], exact_inner_indices[i]);
    }
    for (std::size_t i{0}; i < outer_indices_size; ++i) {
        CHECK_EQ(outer_indices[i], exact_outer_indices[i]);
    }
}

TEST_CASE("Resize") {
    CscMatrix csc_matrix(5, 5);

    csc_matrix.updateCoeff(4, 4, 8.0);
    csc_matrix.updateCoeff(1, 4, 17.0);
    csc_matrix.updateCoeff(2, 3, 1.0);
    csc_matrix.updateCoeff(4, 2, 14.0);
    csc_matrix.updateCoeff(2, 1, 5.0);
    csc_matrix.updateCoeff(0, 1, 3.0);
    csc_matrix.updateCoeff(2, 0, 7.0);
    csc_matrix.updateCoeff(1, 0, 22.0);
    csc_matrix.updateCoeff(5, 5, 1.0);

    SUBCASE("EnlargeSize") {
        csc_matrix.resize(10, 10);

        const std::size_t values_size{csc_matrix.values().size()};
        const std::size_t inner_indices_size{csc_matrix.innerIndices().size()};
        const std::size_t outer_indices_size{csc_matrix.outerIndices().size()};
        const std::vector<double>& values{csc_matrix.values()};
        const std::vector<int>& inner_indices{csc_matrix.innerIndices()};
        const std::vector<int>& outer_indices{csc_matrix.outerIndices()};

        const std::vector<double> exact_values{22.0, 7.0, 3.0, 5.0, 14.0, 1.0, 17.0, 8.0};
        const std::vector<int> exact_inner_indices{1, 2, 0, 2, 4, 2, 1, 4};
        const std::vector<int> exact_outer_indices{0, 2, 4, 5, 6, 8, 8, 8, 8, 8, 8};

        CHECK_EQ(values_size, exact_values.size());
        CHECK_EQ(inner_indices_size, exact_inner_indices.size());
        CHECK_EQ(outer_indices_size, exact_outer_indices.size());

        for (std::size_t i{0}; i < values_size; ++i) {
            CHECK_EQ(values[i], exact_values[i]);
        }
        for (std::size_t i{0}; i < inner_indices_size; ++i) {
            CHECK_EQ(inner_indices[i], exact_inner_indices[i]);
        }
        for (std::size_t i{0}; i < outer_indices_size; ++i) {
            CHECK_EQ(outer_indices[i], exact_outer_indices[i]);
        }
    }

    SUBCASE("DescendSize") {
        csc_matrix.resize(3, 3);

        const std::size_t values_size{csc_matrix.values().size()};
        const std::size_t inner_indices_size{csc_matrix.innerIndices().size()};
        const std::size_t outer_indices_size{csc_matrix.outerIndices().size()};
        const std::vector<double>& values{csc_matrix.values()};
        const std::vector<int>& inner_indices{csc_matrix.innerIndices()};
        const std::vector<int>& outer_indices{csc_matrix.outerIndices()};

        const std::vector<double> exact_values{22.0, 7.0, 3.0, 5.0};
        const std::vector<int> exact_inner_indices{1, 2, 0, 2};
        const std::vector<int> exact_outer_indices{0, 2, 4, 4};

        CHECK_EQ(values_size, exact_values.size());
        CHECK_EQ(inner_indices_size, exact_inner_indices.size());
        CHECK_EQ(outer_indices_size, exact_outer_indices.size());

        for (std::size_t i{0}; i < values_size; ++i) {
            CHECK_EQ(values[i], exact_values[i]);
        }
        for (std::size_t i{0}; i < inner_indices_size; ++i) {
            CHECK_EQ(inner_indices[i], exact_inner_indices[i]);
        }
        for (std::size_t i{0}; i < outer_indices_size; ++i) {
            CHECK_EQ(outer_indices[i], exact_outer_indices[i]);
        }
    }
}

TEST_CASE("Serialization") {
    CscMatrix csc_matrix(5, 5);

    csc_matrix.updateCoeff(4, 4, 8.0);
    csc_matrix.updateCoeff(1, 4, 17.0);
    csc_matrix.updateCoeff(2, 3, 1.0);
    csc_matrix.updateCoeff(4, 2, 14.0);
    csc_matrix.updateCoeff(2, 1, 5.0);
    csc_matrix.updateCoeff(0, 1, 3.0);
    csc_matrix.updateCoeff(2, 0, 7.0);
    csc_matrix.updateCoeff(1, 0, 22.0);
    csc_matrix.updateCoeff(5, 5, 1.0);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << csc_matrix;

    CscMatrix other_csc_matrix;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_csc_matrix;

    const std::size_t values_size{other_csc_matrix.values().size()};
    const std::size_t inner_indices_size{other_csc_matrix.innerIndices().size()};
    const std::size_t outer_indices_size{other_csc_matrix.outerIndices().size()};
    const std::vector<double>& values{other_csc_matrix.values()};
    const std::vector<int>& inner_indices{other_csc_matrix.innerIndices()};
    const std::vector<int>& outer_indices{other_csc_matrix.outerIndices()};

    const std::vector<double> exact_values{22.0, 7.0, 3.0, 5.0, 14.0, 1.0, 17.0, 8.0};
    const std::vector<int> exact_inner_indices{1, 2, 0, 2, 4, 2, 1, 4};
    const std::vector<int> exact_outer_indices{0, 2, 4, 5, 6, 8};

    CHECK_EQ(values_size, exact_values.size());
    CHECK_EQ(inner_indices_size, exact_inner_indices.size());
    CHECK_EQ(outer_indices_size, exact_outer_indices.size());

    for (std::size_t i{0}; i < values_size; ++i) {
        CHECK_EQ(values[i], exact_values[i]);
    }
    for (std::size_t i{0}; i < inner_indices_size; ++i) {
        CHECK_EQ(inner_indices[i], exact_inner_indices[i]);
    }
    for (std::size_t i{0}; i < outer_indices_size; ++i) {
        CHECK_EQ(outer_indices[i], exact_outer_indices[i]);
    }
}

} // namespace math
} // namespace tiny_pnc
