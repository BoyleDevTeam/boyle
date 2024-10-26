/**
 * @file eigen_serialization.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-27
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <vector>

#include "boost/serialization/array.hpp"
#include "boost/serialization/split_free.hpp"
#include "boost/serialization/vector.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"

namespace boost::serialization {

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline auto save(
    auto& archive, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    const int rows = m.rows();
    const int cols = m.cols();
    archive & rows;
    archive & cols;
    archive& boost::serialization::make_array(m.data(), rows * cols);
    return;
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline auto load(
    auto& archive, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    int rows, cols;
    archive & rows;
    archive & cols;
    m.resize(rows, cols);
    archive& boost::serialization::make_array(m.data(), rows * cols);
    return;
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    split_free(archive, m, version);
    return;
}

template <typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline auto save(
    auto& archive, const Eigen::Triplet<_Scalar, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & m.row();
    archive & m.col();
    archive & m.value();
    return;
}

template <typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline auto load(
    auto& archive, Eigen::Triplet<_Scalar, _Index>& m, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    int row, col;
    _Scalar value;
    archive & row;
    archive & col;
    archive & value;
    m = Eigen::Triplet<_Scalar, _Index>(row, col, value);
    return;
}

template <typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, Eigen::Triplet<_Scalar, _Index>& m, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    split_free(archive, m, version);
    return;
}

template <typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline auto save(
    auto& archive, const Eigen::SparseMatrix<_Scalar, _Options, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    const int innerSize = m.innerSize();
    const int outerSize = m.outerSize();
    std::vector<Eigen::Triplet<_Scalar>> triplets;
    triplets.reserve(outerSize);
    for (int i{0}; i < outerSize; ++i) {
        for (typename Eigen::SparseMatrix<_Scalar, _Options, _Index>::InnerIterator it(m, i); it;
             ++it) {
            triplets.push_back(Eigen::Triplet<_Scalar>(it.row(), it.col(), it.value()));
        }
    }
    archive & innerSize;
    archive & outerSize;
    archive & triplets;
    return;
}

template <typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline auto load(
    auto& archive, Eigen::SparseMatrix<_Scalar, _Options, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    int innerSize;
    int outerSize;
    archive & innerSize;
    archive & outerSize;
    const int rows = m.IsRowMajor ? outerSize : innerSize;
    const int cols = m.IsRowMajor ? innerSize : outerSize;
    m.resize(rows, cols);
    std::vector<Eigen::Triplet<_Scalar>> triplets;
    archive & triplets;
    m.setFromTriplets(triplets.begin(), triplets.end());
    return;
}

template <typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, Eigen::SparseMatrix<_Scalar, _Options, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    split_free(archive, m, version);
    return;
}

} // namespace boost::serialization
