/**
 * @file lil_matrix_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-28
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/sparse/lil_matrix.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "boyle/math/sparse/coo_matrix.hpp"
#include "boyle/math/sparse/csc_matrix.hpp"
#include "boyle/math/sparse/csr_matrix.hpp"
#include "boyle/math/sparse/dok_matrix.hpp"
#include "boyle/math/sparse/dok_vector.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Basic") {
    LilMatrix<double, int> lil_matrix(4, 8);

    DokVector<double, int> row_vector_0{{3, 1.5}, {1, 345.2}, {0, 567.4}, {5, 0.0}};
    DokVector<double, int> row_vector_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, row_vector_0);
    lil_matrix.updateRow(3, row_vector_3);

    CHECK_EQ(lil_matrix.nrows(), 4);
    CHECK_EQ(lil_matrix.ncols(), 8);
    CHECK_EQ(lil_matrix.nnzs(), 8);

    CHECK_EQ(lil_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(lil_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(lil_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(lil_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(lil_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 2), 1.57);
    CHECK_EQ(lil_matrix.coeff(3, 4), -35.2);
    CHECK_EQ(lil_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 3), 3775.0);

    CHECK_EQ(lil_matrix[0, 3], 1.5);
    CHECK_EQ(lil_matrix[0, 1], 345.2);
    CHECK_EQ(lil_matrix[0, 0], 567.4);
    CHECK_EQ(lil_matrix[0, 5], 0.0);
    CHECK_EQ(lil_matrix[0, 2], 0.0);
    CHECK_EQ(lil_matrix[3, 2], 1.57);
    CHECK_EQ(lil_matrix[3, 4], -35.2);
    CHECK_EQ(lil_matrix[3, 1], 0.0);
    CHECK_EQ(lil_matrix[3, 3], 3775.0);
}

TEST_CASE("Resize") {
    LilMatrix<double, int> lil_matrix(4, 8);

    DokVector<double, int> row_vector_0{{3, 1.5}, {1, 345.2}, {0, 567.4}, {5, 0.0}};
    DokVector<double, int> row_vector_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, row_vector_0);
    lil_matrix.updateRow(3, row_vector_3);

    lil_matrix.resize(2, 4);

    CHECK_EQ(lil_matrix.nrows(), 2);
    CHECK_EQ(lil_matrix.ncols(), 4);
    CHECK_EQ(lil_matrix.nnzs(), 3);

    CHECK_EQ(lil_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(lil_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(lil_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(lil_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(lil_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 2), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 4), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 3), 0.0);

    CHECK_EQ(lil_matrix[0, 3], 1.5);
    CHECK_EQ(lil_matrix[0, 1], 345.2);
    CHECK_EQ(lil_matrix[0, 0], 567.4);
    CHECK_EQ(lil_matrix[0, 5], 0.0);
    CHECK_EQ(lil_matrix[0, 2], 0.0);
    CHECK_EQ(lil_matrix[3, 2], 0.0);
    CHECK_EQ(lil_matrix[3, 4], 0.0);
    CHECK_EQ(lil_matrix[3, 1], 0.0);
    CHECK_EQ(lil_matrix[3, 3], 0.0);
}

TEST_SUITE("Conversion") {
    TEST_CASE_TEMPLATE(
        "FromLilMatrix", T, DokMatrix<double, int>, CooMatrix<double, int>, CscMatrix<double, int>,
        CsrMatrix<double, int>
    ) {
        LilMatrix<double, int> lil_matrix(4, 8);

        lil_matrix.updateCoeff(0, 3, 1.5);
        lil_matrix.updateCoeff(0, 1, 345.2);
        lil_matrix.updateCoeff(0, 0, 567.4);
        lil_matrix.updateCoeff(0, 5, 0.0);
        lil_matrix.updateCoeff(3, 2, 1.57);
        lil_matrix.updateCoeff(3, 4, -35.2);
        lil_matrix.updateCoeff(3, 1, 0.0);
        lil_matrix.updateCoeff(3, 3, 3775.0);

        T sparse_matrix{lil_matrix};

        CHECK_EQ(sparse_matrix.nrows(), 4);
        CHECK_EQ(sparse_matrix.ncols(), 8);
        CHECK_EQ(sparse_matrix.nnzs(), 8);

        CHECK_EQ(sparse_matrix.coeff(0, 3), 1.5);
        CHECK_EQ(sparse_matrix.coeff(0, 1), 345.2);
        CHECK_EQ(sparse_matrix.coeff(0, 0), 567.4);
        CHECK_EQ(sparse_matrix.coeff(0, 5), 0.0);
        CHECK_EQ(sparse_matrix.coeff(0, 2), 0.0);
        CHECK_EQ(sparse_matrix.coeff(3, 2), 1.57);
        CHECK_EQ(sparse_matrix.coeff(3, 4), -35.2);
        CHECK_EQ(sparse_matrix.coeff(3, 1), 0.0);
        CHECK_EQ(sparse_matrix.coeff(3, 3), 3775.0);

        CHECK_EQ(sparse_matrix[0, 3], 1.5);
        CHECK_EQ(sparse_matrix[0, 1], 345.2);
        CHECK_EQ(sparse_matrix[0, 0], 567.4);
        CHECK_EQ(sparse_matrix[0, 5], 0.0);
        CHECK_EQ(sparse_matrix[0, 2], 0.0);
        CHECK_EQ(sparse_matrix[3, 2], 1.57);
        CHECK_EQ(sparse_matrix[3, 4], -35.2);
        CHECK_EQ(sparse_matrix[3, 1], 0.0);
        CHECK_EQ(sparse_matrix[3, 3], 3775.0);
    }

    TEST_CASE_TEMPLATE(
        "ToLilMatrix", T, DokMatrix<double, int>, CooMatrix<double, int>, CscMatrix<double, int>,
        CsrMatrix<double, int>
    ) {
        T sparse_matrix(4, 8);

        sparse_matrix.updateCoeff(0, 3, 1.5);
        sparse_matrix.updateCoeff(0, 1, 345.2);
        sparse_matrix.updateCoeff(0, 0, 567.4);
        sparse_matrix.updateCoeff(0, 5, 0.0);
        sparse_matrix.updateCoeff(3, 2, 1.57);
        sparse_matrix.updateCoeff(3, 4, -35.2);
        sparse_matrix.updateCoeff(3, 1, 0.0);
        sparse_matrix.updateCoeff(3, 3, 3775.0);

        LilMatrix<double, int> lil_matrix{sparse_matrix};

        CHECK_EQ(lil_matrix.nrows(), 4);
        CHECK_EQ(lil_matrix.ncols(), 8);
        CHECK_EQ(lil_matrix.nnzs(), 8);

        CHECK_EQ(lil_matrix.coeff(0, 3), 1.5);
        CHECK_EQ(lil_matrix.coeff(0, 1), 345.2);
        CHECK_EQ(lil_matrix.coeff(0, 0), 567.4);
        CHECK_EQ(lil_matrix.coeff(0, 5), 0.0);
        CHECK_EQ(lil_matrix.coeff(0, 2), 0.0);
        CHECK_EQ(lil_matrix.coeff(3, 2), 1.57);
        CHECK_EQ(lil_matrix.coeff(3, 4), -35.2);
        CHECK_EQ(lil_matrix.coeff(3, 1), 0.0);
        CHECK_EQ(lil_matrix.coeff(3, 3), 3775.0);

        CHECK_EQ(lil_matrix[0, 3], 1.5);
        CHECK_EQ(lil_matrix[0, 1], 345.2);
        CHECK_EQ(lil_matrix[0, 0], 567.4);
        CHECK_EQ(lil_matrix[0, 5], 0.0);
        CHECK_EQ(lil_matrix[0, 2], 0.0);
        CHECK_EQ(lil_matrix[3, 2], 1.57);
        CHECK_EQ(lil_matrix[3, 4], -35.2);
        CHECK_EQ(lil_matrix[3, 1], 0.0);
        CHECK_EQ(lil_matrix[3, 3], 3775.0);
    }
}

TEST_CASE("Serialization") {
    LilMatrix<double, int> lil_matrix(4, 8);

    DokVector<double, int> row_vector_0{{3, 1.5}, {1, 345.2}, {0, 567.4}, {5, 0.0}};
    DokVector<double, int> row_vector_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, std::move(row_vector_0));
    lil_matrix.updateRow(3, row_vector_3);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << lil_matrix;

    LilMatrix<double, int> other_lil_matrix;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_lil_matrix;

    CHECK_EQ(other_lil_matrix.nrows(), 4);
    CHECK_EQ(other_lil_matrix.ncols(), 8);
    CHECK_EQ(other_lil_matrix.nnzs(), 8);

    CHECK_EQ(other_lil_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(other_lil_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(other_lil_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(other_lil_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(3, 2), 1.57);
    CHECK_EQ(other_lil_matrix.coeff(3, 4), -35.2);
    CHECK_EQ(other_lil_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(3, 3), 3775.0);

    CHECK_EQ(other_lil_matrix[0, 3], 1.5);
    CHECK_EQ(other_lil_matrix[0, 1], 345.2);
    CHECK_EQ(other_lil_matrix[0, 0], 567.4);
    CHECK_EQ(other_lil_matrix[0, 5], 0.0);
    CHECK_EQ(other_lil_matrix[0, 2], 0.0);
    CHECK_EQ(other_lil_matrix[3, 2], 1.57);
    CHECK_EQ(other_lil_matrix[3, 4], -35.2);
    CHECK_EQ(other_lil_matrix[3, 1], 0.0);
    CHECK_EQ(other_lil_matrix[3, 3], 3775.0);
}

} // namespace boyle::math
