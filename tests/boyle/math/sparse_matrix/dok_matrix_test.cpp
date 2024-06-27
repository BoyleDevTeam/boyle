/**
 * @file dok_matrix_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-25
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/sparse_matrix/dok_matrix.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "boyle/math/sparse_matrix/coo_matrix.hpp"
#include "boyle/math/sparse_matrix/csc_matrix.hpp"
#include "boyle/math/sparse_matrix/csr_matrix.hpp"
#include "boyle/math/sparse_matrix/lil_matrix.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Basic") {
    DokMatrix dok_matrix(8, 8);

    CHECK_EQ(dok_matrix.nrows(), 8);
    CHECK_EQ(dok_matrix.ncols(), 8);

    dok_matrix.updateCoeff(2, 5, 5.0);
    dok_matrix.updateCoeff(7, 3, 2874.0843);
    dok_matrix.updateCoeff(0, 6, -408.876);
    dok_matrix.updateCoeff(4, 6, 0.0);
    dok_matrix.updateCoeff(8, 0, 1.5);

    CHECK_EQ(dok_matrix.coeff(2, 5), 5.0);
    CHECK_EQ(dok_matrix.coeff(7, 3), 2874.0843);
    CHECK_EQ(dok_matrix.coeff(0, 6), -408.876);
    CHECK_EQ(dok_matrix.coeff(4, 6), 0.0);
    CHECK_EQ(dok_matrix.coeff(8, 0), 0.0);

    CHECK_EQ(dok_matrix.coeff(1, 3), 0.0);
    CHECK_EQ(dok_matrix.coeff(2, 4), 0.0);

    std::size_t nnzs = dok_matrix.nnzs();

    CHECK_EQ(nnzs, 3);

    dok_matrix.updateCoeff(7, 3, 28234.0843);
    dok_matrix.updateCoeff(2, 5, 7.0);
    dok_matrix.updateCoeff(0, 6, 408.876);
    dok_matrix.updateCoeff(4, 6, 1.0);

    nnzs = dok_matrix.nnzs();
    CHECK_EQ(nnzs, 4);

    CHECK_EQ(dok_matrix.coeff(2, 5), 7.0);
    CHECK_EQ(dok_matrix.coeff(7, 3), 28234.0843);
    CHECK_EQ(dok_matrix.coeff(0, 6), 408.876);
    CHECK_EQ(dok_matrix.coeff(4, 6), 1.0);

    dok_matrix.clear();

    nnzs = dok_matrix.nnzs();
    CHECK_EQ(nnzs, 0);
}

TEST_CASE("Resize") {
    DokMatrix dok_matrix(8, 8);

    CHECK_EQ(dok_matrix.nrows(), 8);
    CHECK_EQ(dok_matrix.ncols(), 8);

    dok_matrix.updateCoeff(2, 5, 5.0);
    dok_matrix.updateCoeff(7, 3, 2874.0843);
    dok_matrix.updateCoeff(0, 6, -408.876);
    dok_matrix.updateCoeff(4, 6, 0.0);
    dok_matrix.updateCoeff(8, 0, 1.5);

    CHECK_EQ(dok_matrix.nnzs(), 3);

    dok_matrix.resize(6, 6);

    CHECK_EQ(dok_matrix.nrows(), 6);
    CHECK_EQ(dok_matrix.ncols(), 6);

    dok_matrix.updateCoeff(2, 5, 5.0);
    dok_matrix.updateCoeff(7, 3, 0.0);
    dok_matrix.updateCoeff(0, 6, 0.0);
    dok_matrix.updateCoeff(4, 6, 0.0);
    dok_matrix.updateCoeff(8, 0, 0.0);

    CHECK_EQ(dok_matrix.nnzs(), 1);
}

TEST_SUITE("Conversion") {
    TEST_CASE_TEMPLATE("FromDokMatrix", T, CooMatrix<double, int>, CscMatrix<double, int>, CsrMatrix<double, int>, LilMatrix<double, int>) {
        DokMatrix<double, int> dok_matrix(4, 8);

        dok_matrix.updateCoeff(0, 3, 1.5);
        dok_matrix.updateCoeff(0, 1, 345.2);
        dok_matrix.updateCoeff(0, 0, 567.4);
        dok_matrix.updateCoeff(0, 5, 0.0);
        dok_matrix.updateCoeff(3, 2, 1.57);
        dok_matrix.updateCoeff(3, 4, -35.2);
        dok_matrix.updateCoeff(3, 1, 0.0);
        dok_matrix.updateCoeff(3, 3, 3775.0);

        T sparse_matrix{dok_matrix};

        CHECK_EQ(sparse_matrix.nrows(), 4);
        CHECK_EQ(sparse_matrix.ncols(), 8);
        CHECK_EQ(sparse_matrix.nnzs(), 6);

        CHECK_EQ(sparse_matrix.coeff(0, 3), 1.5);
        CHECK_EQ(sparse_matrix.coeff(0, 1), 345.2);
        CHECK_EQ(sparse_matrix.coeff(0, 0), 567.4);
        CHECK_EQ(sparse_matrix.coeff(0, 5), 0.0);
        CHECK_EQ(sparse_matrix.coeff(0, 2), 0.0);
        CHECK_EQ(sparse_matrix.coeff(3, 2), 1.57);
        CHECK_EQ(sparse_matrix.coeff(3, 4), -35.2);
        CHECK_EQ(sparse_matrix.coeff(3, 1), 0.0);
        CHECK_EQ(sparse_matrix.coeff(3, 3), 3775.0);
    }

    TEST_CASE_TEMPLATE("ToDokMatrix", T, CooMatrix<double, int>, CscMatrix<double, int>, CsrMatrix<double, int>, LilMatrix<double, int>) {
        T sparse_matrix(4, 8);

        sparse_matrix.updateCoeff(0, 3, 1.5);
        sparse_matrix.updateCoeff(0, 1, 345.2);
        sparse_matrix.updateCoeff(0, 0, 567.4);
        sparse_matrix.updateCoeff(0, 5, 0.0);
        sparse_matrix.updateCoeff(3, 2, 1.57);
        sparse_matrix.updateCoeff(3, 4, -35.2);
        sparse_matrix.updateCoeff(3, 1, 0.0);
        sparse_matrix.updateCoeff(3, 3, 3775.0);

        DokMatrix<double, int> dok_matrix{sparse_matrix};

        CHECK_EQ(dok_matrix.nrows(), 4);
        CHECK_EQ(dok_matrix.ncols(), 8);
        CHECK_EQ(dok_matrix.nnzs(), 6);

        CHECK_EQ(dok_matrix.coeff(0, 3), 1.5);
        CHECK_EQ(dok_matrix.coeff(0, 1), 345.2);
        CHECK_EQ(dok_matrix.coeff(0, 0), 567.4);
        CHECK_EQ(dok_matrix.coeff(0, 5), 0.0);
        CHECK_EQ(dok_matrix.coeff(0, 2), 0.0);
        CHECK_EQ(dok_matrix.coeff(3, 2), 1.57);
        CHECK_EQ(dok_matrix.coeff(3, 4), -35.2);
        CHECK_EQ(dok_matrix.coeff(3, 1), 0.0);
        CHECK_EQ(dok_matrix.coeff(3, 3), 3775.0);
    }
}

TEST_CASE("Serialization") {
    DokMatrix dok_matrix(8, 8);

    dok_matrix.updateCoeff(2, 5, 5.0);
    dok_matrix.updateCoeff(7, 3, 2874.0843);
    dok_matrix.updateCoeff(0, 6, -408.876);
    dok_matrix.updateCoeff(4, 6, 0.0);
    dok_matrix.updateCoeff(8, 0, 1.5);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << dok_matrix;

    DokMatrix other_dok_matrix;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_dok_matrix;

    CHECK_EQ(other_dok_matrix.nrows(), 8);
    CHECK_EQ(other_dok_matrix.ncols(), 8);
    CHECK_EQ(other_dok_matrix.nnzs(), 3);

    CHECK_EQ(other_dok_matrix.coeff(2, 5), 5.0);
    CHECK_EQ(other_dok_matrix.coeff(7, 3), 2874.0843);
    CHECK_EQ(other_dok_matrix.coeff(0, 6), -408.876);
    CHECK_EQ(other_dok_matrix.coeff(4, 6), 0.0);
    CHECK_EQ(other_dok_matrix.coeff(8, 0), 0.0);
}

} // namespace boyle::math
