/**
 * @file dense_degenerate_trait.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-01-18
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"

namespace boyle::math::detail {

template <typename T, ::boyle::math::MatrixOrder Order = ::boyle::math::MatrixOrder::COL_MAJOR>
struct DenseDegenerateTrait final {
    static_assert(false, "MatrixDegenerationTrait is not implemented for this type.");
};

template <typename T, ::boyle::math::MatrixOrder Order = ::boyle::math::MatrixOrder::COL_MAJOR>
using DenseDegenerateTraitT = typename DenseDegenerateTrait<T, Order>::type;

template <ScalarArithmetic Scalar>
struct DenseDegenerateTrait<Scalar, ::boyle::math::MatrixOrder::COL_MAJOR> final {
    using type = Scalar;
};

template <ScalarArithmetic Scalar>
struct DenseDegenerateTrait<Scalar, ::boyle::math::MatrixOrder::ROW_MAJOR> final {
    using type = Scalar;
};

template <
    ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols, ::boyle::math::MatrixOrder Order>
struct DenseDegenerateTrait<
    Matrix<Scalar, NRows, NCols, Order>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = Vector<Scalar, NRows>;
};

template <
    ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols, ::boyle::math::MatrixOrder Order>
struct DenseDegenerateTrait<
    Matrix<Scalar, NRows, NCols, Order>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = Vector<Scalar, NCols>;
};

template <ScalarArithmetic Scalar, ::boyle::math::MatrixOrder Order, Allocatory Alloc>
struct DenseDegenerateTrait<MatrixX<Scalar, Order, Alloc>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = VectorX<Scalar, Alloc>;
};

template <ScalarArithmetic Scalar, ::boyle::math::MatrixOrder Order, Allocatory Alloc>
struct DenseDegenerateTrait<MatrixX<Scalar, Order, Alloc>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = VectorX<Scalar, Alloc>;
};

template <ScalarArithmetic Scalar, std::size_t N>
struct DenseDegenerateTrait<Vector<Scalar, N>, ::boyle::math::MatrixOrder::COL_MAJOR> final {
    using type = Scalar;
};

template <ScalarArithmetic Scalar, std::size_t N>
struct DenseDegenerateTrait<Vector<Scalar, N>, ::boyle::math::MatrixOrder::ROW_MAJOR> final {
    using type = Scalar;
};

template <ScalarArithmetic Scalar, Allocatory Alloc>
struct DenseDegenerateTrait<VectorX<Scalar, Alloc>, ::boyle::math::MatrixOrder::COL_MAJOR> final {
    using type = Scalar;
};

template <ScalarArithmetic Scalar, Allocatory Alloc>
struct DenseDegenerateTrait<VectorX<Scalar, Alloc>, ::boyle::math::MatrixOrder::ROW_MAJOR> final {
    using type = Scalar;
};

} // namespace boyle::math::detail
