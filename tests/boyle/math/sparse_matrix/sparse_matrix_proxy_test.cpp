/**
 * @file sparse_matrix_proxy_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/sparse_matrix/sparse_matrix_proxy.hpp"

#include <vector>

#include "boyle/math/sparse_matrix/coo_matrix.hpp"
#include "boyle/math/sparse_matrix/csc_matrix.hpp"
#include "boyle/math/sparse_matrix/csr_matrix.hpp"
#include "boyle/math/sparse_matrix/dok_matrix.hpp"
#include "boyle/math/sparse_matrix/lil_matrix.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

TEST_CASE("Polymorphism") {
    std::vector<SparseMatrixProxy<double, int>> matrices{};

    DokMatrix dok_matrix(4, 8);
    dok_matrix.updateCoeff(0, 3, 1.5);
    dok_matrix.updateCoeff(0, 1, 345.2);
    dok_matrix.updateCoeff(0, 0, 567.4);
    dok_matrix.updateCoeff(0, 5, 0.0);
    dok_matrix.updateCoeff(3, 2, 1.57);
    dok_matrix.updateCoeff(3, 4, -35.2);
    dok_matrix.updateCoeff(3, 1, 0.0);
    dok_matrix.updateCoeff(3, 3, 3775.0);

    LilMatrix lil_matrix(4, 8);
    lil_matrix.updateCoeff(0, 3, 1.5);
    lil_matrix.updateCoeff(0, 1, 345.2);
    lil_matrix.updateCoeff(0, 0, 567.4);
    lil_matrix.updateCoeff(0, 5, 0.0);
    lil_matrix.updateCoeff(3, 2, 1.57);
    lil_matrix.updateCoeff(3, 4, -35.2);
    lil_matrix.updateCoeff(3, 1, 0.0);
    lil_matrix.updateCoeff(3, 3, 3775.0);

    CooMatrix coo_matrix(4, 8);
    coo_matrix.updateCoeff(0, 3, 1.5);
    coo_matrix.updateCoeff(0, 1, 345.2);
    coo_matrix.updateCoeff(0, 0, 567.4);
    coo_matrix.updateCoeff(0, 5, 0.0);
    coo_matrix.updateCoeff(3, 2, 1.57);
    coo_matrix.updateCoeff(3, 4, -35.2);
    coo_matrix.updateCoeff(3, 1, 0.0);
    coo_matrix.updateCoeff(3, 3, 3775.0);

    CscMatrix csc_matrix(4, 8);
    csc_matrix.updateCoeff(0, 3, 1.5);
    csc_matrix.updateCoeff(0, 1, 345.2);
    csc_matrix.updateCoeff(0, 0, 567.4);
    csc_matrix.updateCoeff(0, 5, 0.0);
    csc_matrix.updateCoeff(3, 2, 1.57);
    csc_matrix.updateCoeff(3, 4, -35.2);
    csc_matrix.updateCoeff(3, 1, 0.0);
    csc_matrix.updateCoeff(3, 3, 3775.0);

    CsrMatrix csr_matrix(4, 8);
    csr_matrix.updateCoeff(0, 3, 1.5);
    csr_matrix.updateCoeff(0, 1, 345.2);
    csr_matrix.updateCoeff(0, 0, 567.4);
    csr_matrix.updateCoeff(0, 5, 0.0);
    csr_matrix.updateCoeff(3, 2, 1.57);
    csr_matrix.updateCoeff(3, 4, -35.2);
    csr_matrix.updateCoeff(3, 1, 0.0);
    csr_matrix.updateCoeff(3, 3, 3775.0);

    matrices.emplace_back(makeSparseMatrixProxy(dok_matrix));
    matrices.emplace_back(makeSparseMatrixProxy(std::move(lil_matrix)));
    matrices.emplace_back(makeSparseMatrixProxy(std::move(coo_matrix)));
    matrices.emplace_back(makeSparseMatrixProxy(std::move(csc_matrix)));
    matrices.emplace_back(makeSparseMatrixProxy(std::move(csr_matrix)));

    for (const SparseMatrixProxy<double, int>& sparse_matrix : matrices) {
        CHECK_EQ(sparse_matrix->nrows(), 4);
        CHECK_EQ(sparse_matrix->ncols(), 8);
        CHECK_EQ(sparse_matrix->nnzs(), 6);

        CHECK_EQ(sparse_matrix->coeff(0, 3), 1.5);
        CHECK_EQ(sparse_matrix->coeff(0, 1), 345.2);
        CHECK_EQ(sparse_matrix->coeff(0, 0), 567.4);
        CHECK_EQ(sparse_matrix->coeff(0, 5), 0.0);
        CHECK_EQ(sparse_matrix->coeff(0, 2), 0.0);
        CHECK_EQ(sparse_matrix->coeff(3, 2), 1.57);
        CHECK_EQ(sparse_matrix->coeff(3, 4), -35.2);
        CHECK_EQ(sparse_matrix->coeff(3, 1), 0.0);
        CHECK_EQ(sparse_matrix->coeff(3, 3), 3775.0);
    }
}

} // namespace boyle::math
