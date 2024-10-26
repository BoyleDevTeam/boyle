/**
 * @file dense_utils.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-04-10
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <iomanip>
#include <ostream>

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"

namespace boyle::math {

template <typename Char, VecArithmetic VectorType>
inline auto operator<<(std::basic_ostream<Char>& os, VectorType&& vector) noexcept
    -> std::basic_ostream<Char>& {
    using value_type = typename std::remove_cvref_t<VectorType>::value_type;
    using size_type = typename std::remove_cvref_t<VectorType>::size_type;
    constexpr size_type kWidth{isComplexArithmeticV<value_type> ? 32 : 16};
    const size_type size{vector.size()};
    os << std::fixed;
    for (size_type i{0}; i < size; ++i) {
        os << std::setw(kWidth) << vector[i] << '\n';
    }
    return os;
}

template <typename Char, MatArithmetic MatrixType>
inline auto operator<<(std::basic_ostream<Char>& os, MatrixType&& matrix) noexcept
    -> std::basic_ostream<Char>& {
    using value_type = typename std::remove_cvref_t<MatrixType>::value_type;
    using size_type = typename std::remove_cvref_t<MatrixType>::size_type;
    constexpr size_type kWidth{isComplexArithmeticV<value_type> ? 32 : 16};
    const size_type nrows{matrix.nrows()}, ncols{matrix.ncols()};
    os << std::fixed;
    for (size_type i{0}; i < nrows; ++i) {
        for (size_type j{0}; j < ncols; ++j) {
            os << std::setw(kWidth) << matrix[i, j];
        }
        os << '\n';
    }
    return os;
}

template <typename Char, ScalarArithmetic Scalar>
inline auto operator<<(std::basic_ostream<Char>& os, VectorView<Scalar> vector_view) noexcept
    -> std::basic_ostream<Char>& {
    using value_type = std::remove_const_t<Scalar>;
    using size_type = typename VectorView<Scalar>::size_type;
    constexpr size_type kWidth{isComplexArithmeticV<value_type> ? 32 : 16};
    const size_type size{vector_view.size()};
    os << std::fixed;
    for (size_type i{0}; i < size; ++i) {
        os << std::setw(kWidth) << vector_view[i] << '\n';
    }
    return os;
}

template <typename Char, ScalarArithmetic Scalar, MatrixOrder Order>
inline auto operator<<(std::basic_ostream<Char>& os, MatrixView<Scalar, Order> matrix_view) noexcept
    -> std::basic_ostream<Char>& {
    using value_type = std::remove_const_t<Scalar>;
    using size_type = typename MatrixView<Scalar, Order>::size_type;
    constexpr size_type kWidth{isComplexArithmeticV<value_type> ? 32 : 16};
    const size_type nrows{matrix_view.nrows()}, ncols{matrix_view.ncols()};
    os << std::fixed;
    for (size_type i{0}; i < nrows; ++i) {
        for (size_type j{0}; j < ncols; ++j) {
            os << std::setw(kWidth) << matrix_view[i, j];
        }
        os << '\n';
    }
    return os;
}

} // namespace boyle::math
