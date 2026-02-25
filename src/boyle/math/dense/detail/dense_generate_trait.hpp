/**
 * @file dense_generate_trait.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-09-07
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include "boyle/math/dense/dense_traits.hpp"

namespace boyle::math::detail {

template <typename T, ::boyle::math::MatrixOrder Order = ::boyle::math::MatrixOrder::COL_MAJOR>
struct DenseGenerateTrait final {
    static_assert(false, "MatrixDegenerationTrait is not implemented for this type.");
};

template <typename T, ::boyle::math::MatrixOrder Order = ::boyle::math::MatrixOrder::COL_MAJOR>
using DenseGenerateTraitT = typename DenseGenerateTrait<T, Order>::type;

template <ScalarArithmetic Scalar, std::size_t N>
struct DenseGenerateTrait<::boyle::math::Vector<Scalar, N>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = ::boyle::math::Matrix<Scalar, N, N, ::boyle::math::MatrixOrder::COL_MAJOR>;
};

template <ScalarArithmetic Scalar, std::size_t N>
struct DenseGenerateTrait<::boyle::math::Vector<Scalar, N>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = ::boyle::math::Matrix<Scalar, N, N, ::boyle::math::MatrixOrder::ROW_MAJOR>;
};

template <ScalarArithmetic Scalar, Allocatory Alloc>
struct DenseGenerateTrait<
    ::boyle::math::VectorX<Scalar, Alloc>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = ::boyle::math::MatrixX<Scalar, ::boyle::math::MatrixOrder::COL_MAJOR, Alloc>;
};

template <ScalarArithmetic Scalar, Allocatory Alloc>
struct DenseGenerateTrait<
    ::boyle::math::VectorX<Scalar, Alloc>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = ::boyle::math::MatrixX<Scalar, ::boyle::math::MatrixOrder::ROW_MAJOR, Alloc>;
};

} // namespace boyle::math::detail
