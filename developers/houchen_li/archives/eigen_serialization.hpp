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

#include "boost/serialization/array_wrapper.hpp"
#include "boost/serialization/split_free.hpp"
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
    int rows;
    int cols;
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

template <typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline auto save(
    auto& archive, const Eigen::SparseMatrix<_Scalar, _Options, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    const _Index rows = m.rows();
    const _Index cols = m.cols();
    const _Index non_zeros = m.nonZeros();
    archive & rows;
    archive & cols;
    archive & non_zeros;
    archive& boost::serialization::make_array(m.outerIndexPtr(), m.outerSize() + 1);
    archive& boost::serialization::make_array(m.innerIndexPtr(), non_zeros);
    archive& boost::serialization::make_array(m.valuePtr(), non_zeros);
    return;
}

template <typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline auto load(
    auto& archive, Eigen::SparseMatrix<_Scalar, _Options, _Index>& m,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    _Index rows;
    _Index cols;
    _Index non_zeros;
    archive & rows;
    archive & cols;
    archive & non_zeros;
    m.resize(rows, cols);
    m.resizeNonZeros(non_zeros);
    archive& boost::serialization::make_array(m.outerIndexPtr(), m.outerSize() + 1);
    archive& boost::serialization::make_array(m.innerIndexPtr(), non_zeros);
    archive& boost::serialization::make_array(m.valuePtr(), non_zeros);
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
