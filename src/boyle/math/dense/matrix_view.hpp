/**
 * @file matrix_view.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-04-01
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <type_traits>

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, MatrixOrder Order>
class MatrixView final {
  public:
    using value_type = typename DenseTraits<MatrixView>::value_type;
    using reference = typename DenseTraits<MatrixView>::reference;
    using const_reference = typename DenseTraits<MatrixView>::const_reference;
    using pointer = typename DenseTraits<MatrixView>::pointer;
    using const_pointer = typename DenseTraits<MatrixView>::const_pointer;
    using size_type = typename DenseTraits<MatrixView>::size_type;
    using difference_type = typename DenseTraits<MatrixView>::difference_type;

    static constexpr MatrixOrder kOrder = Order;

    constexpr MatrixView() noexcept = delete;
    constexpr MatrixView(const MatrixView& other) noexcept = default;
    constexpr auto operator=(const MatrixView& other) noexcept -> MatrixView& = default;
    constexpr MatrixView(MatrixView&& other) noexcept = default;
    constexpr auto operator=(MatrixView&&) noexcept -> MatrixView& = default;
    constexpr ~MatrixView() noexcept = default;

    [[using gnu: always_inline]]
    constexpr explicit MatrixView(pointer data, size_type nrows, size_type ncols) noexcept
        : m_data{data}, m_nrows{nrows}, m_ncols{ncols},
          m_stride{kOrder == MatrixOrder::COL_MAJOR ? nrows : ncols} {}

    [[using gnu: always_inline]]
    constexpr explicit MatrixView(
        pointer data, size_type rows, size_type cols, size_type stride
    ) noexcept
        : m_data{data}, m_nrows{rows}, m_ncols{cols}, m_stride{stride} {}

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    constexpr auto nrows() const noexcept -> size_type {
        return m_nrows;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    constexpr auto ncols() const noexcept -> size_type {
        return m_ncols;
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    static constexpr auto order() noexcept -> MatrixOrder {
        return kOrder;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    constexpr auto size() const noexcept -> size_type {
        return m_nrows * m_ncols;
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    constexpr auto stride() const noexcept -> size_type {
        return m_stride;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto data() noexcept -> pointer {
        return m_data;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto data() const noexcept -> const_pointer {
        return m_data;
    }

    [[using gnu: always_inline, hot]]
    constexpr auto fill(const_reference value) noexcept -> void {
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                operator[](i, j) = value;
            }
        }
        return;
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    constexpr auto empty() noexcept -> bool {
        return m_nrows * m_ncols == 0;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type row, size_type col) noexcept -> reference {
        return m_data[offset(row, col)];
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type row, size_type col) const noexcept -> const_reference {
        return m_data[offset(row, col)];
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type i) noexcept -> reference {
        return m_data[offset(i)];
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type i) const noexcept -> const_reference {
        return m_data[offset(i)];
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    constexpr auto coeff(size_type row, size_type col) const noexcept(!BOYLE_CHECK_PARAMS)
        -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (row >= m_nrows || col >= m_ncols) {
            throw std::out_of_range("Matrix index out of range");
        }
#endif
        return operator[](row, col);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto updateCoeff(size_type row, size_type col, const_reference value) noexcept(
        !BOYLE_CHECK_PARAMS
    ) -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (row >= m_nrows || col >= m_ncols) {
            throw std::out_of_range("Matrix index out of range");
        }
#endif
        operator[](row, col) = value;
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    constexpr auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Matrix index out of range.");
        }
#endif
        return operator[](i);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto updateCoeff(size_type i, const_reference value) const noexcept(
        !BOYLE_CHECK_PARAMS
    ) -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Matrix index out of range.");
        }
#endif
        operator[](i) = value;
        return;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view(
        std::pair<size_type, size_type> row_range, std::pair<size_type, size_type> col_range
    ) noexcept -> MatrixView<value_type, kOrder> {
        row_range.second = std::min(row_range.second, m_nrows);
        col_range.second = std::min(col_range.second, m_ncols);
        return MatrixView<value_type, kOrder>(
            data() + offset(row_range.first, col_range.first), row_range.second - row_range.first,
            col_range.second - col_range.first, stride()
        );
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view(
        std::pair<size_type, size_type> row_range, std::pair<size_type, size_type> col_range
    ) const noexcept -> MatrixView<const value_type, kOrder> {
        row_range.second = std::min(row_range.second, m_nrows);
        col_range.second = std::min(col_range.second, m_ncols);
        return MatrixView<const value_type, kOrder>(
            data() + offset(row_range.first, col_range.first), row_range.second - row_range.first,
            col_range.second - col_range.first, stride()

        );
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view(std::pair<size_type, size_type> range) noexcept(!BOYLE_CHECK_PARAMS)
        -> MatrixView<value_type, kOrder> {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != 1 && m_ncols != 1) {
            throw std::invalid_argument("vector does not have setIdentity method.");
        }
#endif
        return view(
            m_ncols == 1 ? range : std::pair<size_type, size_type>{0, 1},
            m_ncols == 1 ? std::pair<size_type, size_type>{0, 1} : range
        );
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view(std::pair<size_type, size_type> range) const noexcept(!BOYLE_CHECK_PARAMS)
        -> MatrixView<const value_type, kOrder> {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != 1 && m_ncols != 1) {
            throw std::invalid_argument("vector does not have setIdentity method.");
        }
#endif
        return view(
            m_ncols == 1 ? range : std::pair<size_type, size_type>{0, 1},
            m_ncols == 1 ? std::pair<size_type, size_type>{0, 1} : range
        );
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view() noexcept -> MatrixView<value_type, kOrder> {
        return view(
            {0, std::numeric_limits<value_type>::max()}, {0, std::numeric_limits<value_type>::max()}
        );
    }

    [[using gnu: pure, always_inline]]
    constexpr auto view() const noexcept -> MatrixView<const value_type, kOrder> {
        return view(
            {0, std::numeric_limits<value_type>::max()}, {0, std::numeric_limits<value_type>::max()}
        );
    }

    [[using gnu: always_inline]]
    constexpr auto setIdentity() noexcept -> void {
        constexpr size_type n{std::min(m_nrows, m_ncols)};
        fill(0.0);
        for (size_type i{0}; i < n; ++i) {
            operator[](i, i) = 1.0;
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator==(const MatrixView& other) const noexcept -> bool {
        if (m_nrows != other.m_nrows || m_ncols != other.m_ncols) [[unlikely]] {
            return false;
        }
        const size_type n{size()};
        for (size_type i{0}; i < n; ++i) {
            if (m_data[i] != other.m_data[i]) {
                return false;
            }
        }
        return true;
    }

    [[using gnu: pure, always_inline, hot]]
    auto identicalTo(MatrixView<value_type, kOrder> obj, double tol = 1E-8) const noexcept -> bool {
        if (m_nrows != obj.nrows() || m_ncols != obj.ncols()) {
            return false;
        }
        const size_type n = size();
        detail::DenseNormTraitT<value_type> test(0.0);
        for (size_type i{0}; i < n; ++i) {
            test = std::max(
                test, std::abs(operator[](i) - obj[i]) /
                          std::max(std::abs(operator[](i)), static_cast<decltype(test)>(1.0))
            );
        }
        return test < tol;
    }

  private:
    [[using gnu: pure, always_inline, leaf, hot]] [[nodiscard]]
    constexpr auto offset(size_type row, size_type col) const noexcept -> size_type {
        if constexpr (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR) {
            return (col * m_stride) + row;
        }
        if constexpr (kOrder == ::boyle::math::MatrixOrder::ROW_MAJOR) {
            return (row * m_stride) + col;
        }
    }

    [[using gnu: pure, always_inline, leaf, hot]] [[nodiscard]]
    constexpr auto offset(size_type i) const noexcept -> size_type {
        if constexpr (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR) {
            return i / m_nrows * m_stride + i % m_nrows;
        }
        if constexpr (kOrder == ::boyle::math::MatrixOrder::ROW_MAJOR) {
            return i / m_ncols * m_stride + i % m_ncols;
        }
    }

    pointer m_data;
    std::size_t m_nrows, m_ncols, m_stride;
};

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
