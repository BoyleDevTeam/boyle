/**
 * @file coo_matrix_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-27
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "math/sparse_matrix/coo_matrix.hpp"

#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "math/sparse_matrix/lil_matrix.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace tiny_pnc {
namespace math {

TEST_CASE("Basic") {
    CooMatrix coo_matrix(8, 8);

    CHECK_EQ(coo_matrix.nrows(), 8);
    CHECK_EQ(coo_matrix.ncols(), 8);

    coo_matrix.updateCoeff(2, 5, 5.0);
    coo_matrix.updateCoeff(7, 3, 2874.0843);
    coo_matrix.updateCoeff(0, 6, -408.876);
    coo_matrix.updateCoeff(4, 6, 0.0);
    coo_matrix.updateCoeff(8, 0, 1.5);

    CHECK_EQ(coo_matrix.coeff(2, 5), 5.0);
    CHECK_EQ(coo_matrix.coeff(7, 3), 2874.0843);
    CHECK_EQ(coo_matrix.coeff(0, 6), -408.876);
    CHECK_EQ(coo_matrix.coeff(4, 6), 0.0);
    CHECK_EQ(coo_matrix.coeff(8, 0), 0.0);

    CHECK_EQ(coo_matrix.coeff(1, 3), 0.0);
    CHECK_EQ(coo_matrix.coeff(2, 4), 0.0);

    std::size_t nnzs = coo_matrix.nnzs();
    const std::vector<double>& values = coo_matrix.values();
    const std::vector<int>& row_indices = coo_matrix.rowIndices();
    const std::vector<int>& col_indices = coo_matrix.colIndices();

    CHECK_EQ(nnzs, 3);

    for (std::size_t i = 0; i < nnzs; ++i) {
        CHECK_EQ(values[i], coo_matrix.coeff(row_indices[i], col_indices[i]));
    }

    coo_matrix.updateCoeff(7, 3, 28234.0843);
    coo_matrix.updateCoeff(2, 5, 7.0);
    coo_matrix.updateCoeff(0, 6, 408.876);
    coo_matrix.updateCoeff(4, 6, 1.0);

    nnzs = coo_matrix.nnzs();
    CHECK_EQ(nnzs, 4);

    CHECK_EQ(coo_matrix.coeff(2, 5), 7.0);
    CHECK_EQ(coo_matrix.coeff(7, 3), 28234.0843);
    CHECK_EQ(coo_matrix.coeff(0, 6), 408.876);
    CHECK_EQ(coo_matrix.coeff(4, 6), 1.0);

    coo_matrix.clear();

    nnzs = coo_matrix.nnzs();
    CHECK_EQ(nnzs, 0);
}

TEST_CASE("Conversion") {
    LilMatrix lil_matrix(4, 8);
    std::vector<double> row_values_0{1.5, 345.2, 567.4, 0.0};
    std::vector<int> row_indices_0{3, 1, 0, 5};

    std::unordered_map<int, double> row_values_map_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, row_values_0, row_indices_0);
    lil_matrix.updateRow(3, row_values_map_3);

    CooMatrix coo_matrix{lil_matrix};

    CHECK_EQ(coo_matrix.nrows(), 4);
    CHECK_EQ(coo_matrix.ncols(), 8);
    CHECK_EQ(coo_matrix.nnzs(), 6);

    CHECK_EQ(coo_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(coo_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(coo_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(coo_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(coo_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(coo_matrix.coeff(3, 2), 1.57);
    CHECK_EQ(coo_matrix.coeff(3, 4), -35.2);
    CHECK_EQ(coo_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(coo_matrix.coeff(3, 3), 3775.0);
}

TEST_CASE("Resize") {
    CooMatrix coo_matrix(8, 8);

    CHECK_EQ(coo_matrix.nrows(), 8);
    CHECK_EQ(coo_matrix.ncols(), 8);

    coo_matrix.updateCoeff(2, 5, 5.0);
    coo_matrix.updateCoeff(7, 3, 2874.0843);
    coo_matrix.updateCoeff(0, 6, -408.876);
    coo_matrix.updateCoeff(4, 6, 0.0);
    coo_matrix.updateCoeff(8, 0, 1.5);

    CHECK_EQ(coo_matrix.nnzs(), 3);

    coo_matrix.resize(6, 6);

    CHECK_EQ(coo_matrix.nrows(), 6);
    CHECK_EQ(coo_matrix.ncols(), 6);

    coo_matrix.updateCoeff(2, 5, 5.0);
    coo_matrix.updateCoeff(7, 3, 0.0);
    coo_matrix.updateCoeff(0, 6, 0.0);
    coo_matrix.updateCoeff(4, 6, 0.0);
    coo_matrix.updateCoeff(8, 0, 0.0);

    CHECK_EQ(coo_matrix.nnzs(), 1);
}

TEST_CASE("Serialization") {
    CooMatrix coo_matrix(8, 8);

    coo_matrix.updateCoeff(2, 5, 5.0);
    coo_matrix.updateCoeff(7, 3, 2874.0843);
    coo_matrix.updateCoeff(0, 6, -408.876);
    coo_matrix.updateCoeff(4, 6, 0.0);
    coo_matrix.updateCoeff(8, 0, 1.5);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << coo_matrix;

    CooMatrix other_coo_matrix;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_coo_matrix;

    CHECK_EQ(other_coo_matrix.nrows(), 8);
    CHECK_EQ(other_coo_matrix.ncols(), 8);
    CHECK_EQ(other_coo_matrix.nnzs(), 3);

    CHECK_EQ(other_coo_matrix.coeff(2, 5), 5.0);
    CHECK_EQ(other_coo_matrix.coeff(7, 3), 2874.0843);
    CHECK_EQ(other_coo_matrix.coeff(0, 6), -408.876);
    CHECK_EQ(other_coo_matrix.coeff(4, 6), 0.0);
    CHECK_EQ(other_coo_matrix.coeff(8, 0), 0.0);
}

} // namespace math
} // namespace tiny_pnc
