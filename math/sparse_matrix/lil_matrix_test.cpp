/**
 * @file lil_matrix_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-28
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "math/sparse_matrix/lil_matrix.hpp"

#include <sstream>
#include <unordered_map>
#include <vector>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#include "math/sparse_matrix/coo_matrix.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace tiny_pnc {
namespace math {

TEST_CASE("Basic") {
    LilMatrix lil_matrix(4, 8);
    std::vector<double> row_values_0{1.5, 345.2, 567.4, 0.0};
    std::vector<int> row_indices_0{3, 1, 0, 5};

    std::unordered_map<int, double> row_values_map_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, row_values_0, row_indices_0);
    lil_matrix.updateRow(3, row_values_map_3);

    CHECK_EQ(lil_matrix.nrows(), 4);
    CHECK_EQ(lil_matrix.ncols(), 8);
    CHECK_EQ(lil_matrix.nnzs(), 6);

    CHECK_EQ(lil_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(lil_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(lil_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(lil_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(lil_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 2), 1.57);
    CHECK_EQ(lil_matrix.coeff(3, 4), -35.2);
    CHECK_EQ(lil_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(lil_matrix.coeff(3, 3), 3775.0);
}

TEST_CASE("Conversion") {
    CooMatrix coo_matrix(8, 8);

    coo_matrix.updateCoeff(2, 5, 5.0);
    coo_matrix.updateCoeff(7, 3, 2874.0843);
    coo_matrix.updateCoeff(0, 6, -408.876);
    coo_matrix.updateCoeff(4, 6, 0.0);
    coo_matrix.updateCoeff(8, 0, 1.5);

    LilMatrix lil_matrix{coo_matrix};

    CHECK_EQ(lil_matrix.nrows(), 8);
    CHECK_EQ(lil_matrix.ncols(), 8);
    CHECK_EQ(lil_matrix.nnzs(), 3);

    CHECK_EQ(lil_matrix.coeff(2, 5), 5.0);
    CHECK_EQ(lil_matrix.coeff(7, 3), 2874.0843);
    CHECK_EQ(lil_matrix.coeff(0, 6), -408.876);
    CHECK_EQ(lil_matrix.coeff(4, 6), 0.0);
    CHECK_EQ(lil_matrix.coeff(8, 0), 0.0);
}

TEST_CASE("Resize") {
    LilMatrix lil_matrix(4, 8);
    std::vector<double> row_values_0{1.5, 345.2, 567.4, 0.0};
    std::vector<int> row_indices_0{3, 1, 0, 5};

    std::unordered_map<int, double> row_values_map_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, row_values_0, row_indices_0);
    lil_matrix.updateRow(3, row_values_map_3);

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
}

TEST_CASE("Serialization") {
    LilMatrix lil_matrix(4, 8);

    std::vector<double> row_values_0{1.5, 345.2, 567.4, 0.0};
    std::vector<int> row_indices_0{3, 1, 0, 5};

    std::unordered_map<int, double> row_values_map_3{{2, 1.57}, {4, -35.2}, {1, 0.0}, {3, 3775.0}};

    lil_matrix.updateRow(0, std::move(row_values_0), std::move(row_indices_0));
    lil_matrix.updateRow(3, row_values_map_3);

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << lil_matrix;

    LilMatrix other_lil_matrix;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_lil_matrix;

    CHECK_EQ(other_lil_matrix.nrows(), 4);
    CHECK_EQ(other_lil_matrix.ncols(), 8);
    CHECK_EQ(other_lil_matrix.nnzs(), 6);

    CHECK_EQ(other_lil_matrix.coeff(0, 3), 1.5);
    CHECK_EQ(other_lil_matrix.coeff(0, 1), 345.2);
    CHECK_EQ(other_lil_matrix.coeff(0, 0), 567.4);
    CHECK_EQ(other_lil_matrix.coeff(0, 5), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(0, 2), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(3, 2), 1.57);
    CHECK_EQ(other_lil_matrix.coeff(3, 4), -35.2);
    CHECK_EQ(other_lil_matrix.coeff(3, 1), 0.0);
    CHECK_EQ(other_lil_matrix.coeff(3, 3), 3775.0);
}

} // namespace math
} // namespace tiny_pnc
