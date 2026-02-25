/**
 * @file dense_partial_pivot_trait.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-01-19
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
struct DensePartialPivotTrait final {
    static_assert(false, "MatrixDegenerationTrait is not implemented for this type.");
};

template <typename T, ::boyle::math::MatrixOrder Order = ::boyle::math::MatrixOrder::COL_MAJOR>
using DensePartialPivotTraitT = typename DensePartialPivotTrait<T, Order>::type;

template <
    ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols, ::boyle::math::MatrixOrder Order>
struct DensePartialPivotTrait<
    Matrix<Scalar, NRows, NCols, Order>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = Vector<std::int32_t, NRows>;
};

template <
    ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols, ::boyle::math::MatrixOrder Order>
struct DensePartialPivotTrait<
    Matrix<Scalar, NRows, NCols, Order>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = Vector<std::int32_t, NCols>;
};

template <ScalarArithmetic Scalar, ::boyle::math::MatrixOrder Order, Allocatory Alloc>
struct DensePartialPivotTrait<MatrixX<Scalar, Order, Alloc>, ::boyle::math::MatrixOrder::COL_MAJOR>
    final {
    using type = VectorX<
        std::int32_t, typename std::allocator_traits<Alloc>::template rebind_alloc<std::int32_t>>;
};

template <ScalarArithmetic Scalar, ::boyle::math::MatrixOrder Order, Allocatory Alloc>
struct DensePartialPivotTrait<MatrixX<Scalar, Order, Alloc>, ::boyle::math::MatrixOrder::ROW_MAJOR>
    final {
    using type = VectorX<
        std::int32_t, typename std::allocator_traits<Alloc>::template rebind_alloc<std::int32_t>>;
};

} // namespace boyle::math::detail
