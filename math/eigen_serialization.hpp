/**
 * @file eigen_serialization.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-27
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include "boost/serialization/array.hpp"
#include "boost/serialization/split_free.hpp"
#include "boost/serialization/vector.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"

namespace boost {
namespace serialization {

template <
    class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline void save(
    Archive& ar, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    const unsigned int version
) noexcept {
    const int rows = m.rows();
    const int cols = m.cols();
    ar& rows;
    ar& cols;
    ar& boost::serialization::make_array(m.data(), rows * cols);
    return;
}

template <
    class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline void load(
    Archive& ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    const unsigned int version
) noexcept {
    int rows, cols;
    ar& rows;
    ar& cols;
    m.resize(rows, cols);
    ar& boost::serialization::make_array(m.data(), rows * cols);
    return;
}

template <
    class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
    const unsigned int version
) noexcept {
    split_free(ar, m, version);
    return;
}

template <class Archive, typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline void save(
    Archive& ar, const Eigen::Triplet<_Scalar, _Index>& m, const unsigned int version
) noexcept {
    ar& m.row();
    ar& m.col();
    ar& m.value();
    return;
}

template <class Archive, typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline void load(
    Archive& ar, Eigen::Triplet<_Scalar, _Index>& m, const unsigned int version
) noexcept {
    int row, col;
    _Scalar value;
    ar& row;
    ar& col;
    ar& value;
    m = Eigen::Triplet<_Scalar, _Index>(row, col, value);
    return;
}

template <class Archive, typename _Scalar, typename _Index>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, Eigen::Triplet<_Scalar, _Index>& m, const unsigned int version
) noexcept {
    split_free(ar, m, version);
    return;
}

template <class Archive, typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline void save(
    Archive& ar, const Eigen::SparseMatrix<_Scalar, _Options, _Index>& m, const unsigned int version
) noexcept {
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
    ar& innerSize;
    ar& outerSize;
    ar& triplets;
    return;
}

template <class Archive, typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline void load(
    Archive& ar, Eigen::SparseMatrix<_Scalar, _Options, _Index>& m, const unsigned int version
) noexcept {
    int innerSize;
    int outerSize;
    ar& innerSize;
    ar& outerSize;
    const int rows = m.IsRowMajor ? outerSize : innerSize;
    const int cols = m.IsRowMajor ? innerSize : outerSize;
    m.resize(rows, cols);
    std::vector<Eigen::Triplet<_Scalar>> triplets;
    ar& triplets;
    m.setFromTriplets(triplets.begin(), triplets.end());
    return;
}

template <class Archive, typename _Scalar, int _Options, typename _Index>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, Eigen::SparseMatrix<_Scalar, _Options, _Index>& m, const unsigned int version
) noexcept {
    split_free(ar, m, version);
    return;
}

} // namespace serialization
} // namespace boost
