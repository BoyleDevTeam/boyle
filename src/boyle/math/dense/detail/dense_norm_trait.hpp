/**
 * @file dense_norm_trait.hpp
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

#include <concepts>

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"

namespace boyle::math::detail {

template <typename T>
struct DenseNormTrait final {
    static_assert(false, "DenseNormTrait is not defined for this type");
};

template <typename T>
using DenseNormTraitT = typename DenseNormTrait<T>::type;

template <std::integral Scalar>
struct DenseNormTrait<Scalar> final {
    using type = Scalar;
};

template <std::floating_point Scalar>
struct DenseNormTrait<Scalar> final {
    using type = Scalar;
};

template <
    std::floating_point Scalar, std::size_t NRows, std::size_t NCols,
    ::boyle::math::MatrixOrder Order>
struct DenseNormTrait<Matrix<Scalar, NRows, NCols, Order>> final {
    using type = Matrix<Scalar, NRows, NCols, Order>;
};

template <std::floating_point Scalar, ::boyle::math::MatrixOrder Order, typename Alloc>
struct DenseNormTrait<MatrixX<Scalar, Order, Alloc>> final {
    using type = MatrixX<Scalar, Order, Alloc>;
};

template <std::floating_point Scalar, std::size_t Size>
struct DenseNormTrait<Vector<Scalar, Size>> final {
    using type = Vector<Scalar, Size>;
};

template <std::floating_point Scalar, typename Alloc>
struct DenseNormTrait<VectorX<Scalar, Alloc>> final {
    using type = VectorX<Scalar, Alloc>;
};

template <ComplexArithmetic Scalar>
struct DenseNormTrait<Scalar> final {
    using type = typename Scalar::value_type;
};

template <
    ComplexArithmetic Scalar, std::size_t NRows, std::size_t NCols,
    ::boyle::math::MatrixOrder Order>
struct DenseNormTrait<Matrix<Scalar, NRows, NCols, Order>> final {
    using type = Matrix<typename Scalar::value_type, NRows, NCols, Order>;
};

template <ComplexArithmetic Scalar, ::boyle::math::MatrixOrder Order, typename Alloc>
struct DenseNormTrait<MatrixX<Scalar, Order, Alloc>> final {
    using type = MatrixX<typename Scalar::value_type, Order, Alloc>;
};

template <ComplexArithmetic Scalar, std::size_t Size>
struct DenseNormTrait<Vector<Scalar, Size>> final {
    using type = Vector<typename Scalar::value_type, Size>;
};

template <ComplexArithmetic Scalar, typename Alloc>
struct DenseNormTrait<VectorX<Scalar, Alloc>> final {
    using type = VectorX<typename Scalar::value_type, Alloc>;
};

} // namespace boyle::math::detail
