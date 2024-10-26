/**
 * @file dense_traits.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-01-08
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <cstdint>
#include <iomanip>
#include <ostream>

#include "boyle/common/utils/aligned_allocator.hpp"
#include "boyle/math/concepts.hpp"

namespace boyle::math {

enum class MatrixOrder : std::int8_t {
    ROW_MAJOR = 101,
    COL_MAJOR = 102
};

template <ScalarArithmetic Scalar, std::size_t N>
class Vector;

template <ScalarArithmetic Scalar, typename Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class VectorX;

template <ScalarArithmetic Scalar>
class VectorView;

template <
    ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols,
    MatrixOrder Order = MatrixOrder::COL_MAJOR>
class Matrix;

template <
    ScalarArithmetic Scalar, MatrixOrder Order = MatrixOrder::COL_MAJOR,
    typename Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
class MatrixX;

template <ScalarArithmetic Scalar, MatrixOrder Order = MatrixOrder::COL_MAJOR>
class MatrixView;

template <typename T>
struct DenseTraits final {
    static_assert(false, "DenseTraits is not implemented for non-Matrix type.");
};

template <ScalarArithmetic Scalar, std::size_t N>
struct DenseTraits<Vector<Scalar, N>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    static constexpr bool is_static = true;
};

template <ScalarArithmetic Scalar, typename Alloc>
struct DenseTraits<VectorX<Scalar, Alloc>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
    static constexpr bool is_static = false;
};

template <ScalarArithmetic Scalar>
struct DenseTraits<VectorView<Scalar>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    static constexpr bool is_static = false;
};

template <ScalarArithmetic Scalar, std::size_t NRows, std::size_t NCols, MatrixOrder Order>
struct DenseTraits<Matrix<Scalar, NRows, NCols, Order>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    static constexpr MatrixOrder kOrder = Order;
    static constexpr bool is_static = true;
};

template <ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
struct DenseTraits<MatrixX<Scalar, Order, Alloc>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;
    static constexpr MatrixOrder kOrder = Order;
    static constexpr bool is_static = false;
};

template <ScalarArithmetic Scalar, MatrixOrder Order>
struct DenseTraits<MatrixView<Scalar, Order>> final {
    using value_type = Scalar;
    using reference = Scalar&;
    using const_reference = const Scalar&;
    using pointer = Scalar*;
    using const_pointer = const Scalar*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    static constexpr MatrixOrder kOrder = Order;
    static constexpr bool is_static = false;
};

} // namespace boyle::math
