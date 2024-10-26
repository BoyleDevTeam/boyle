/**
 * @file matrixx.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-04-02
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <complex>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef BOYLE_USE_BLAS_LAPACK
#include "cblas.h"
#endif

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/dense/matrix_view.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
class MatrixX final {
  public:
    using value_type = Scalar;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = Alloc;

    static constexpr MatrixOrder kOrder = Order;

    MatrixX() noexcept = default;
    MatrixX(const MatrixX& other) = default;
    auto operator=(const MatrixX& other) -> MatrixX& = default;
    MatrixX(MatrixX&& other) noexcept = default;
    auto operator=(MatrixX&& other) noexcept -> MatrixX& = default;
    ~MatrixX() noexcept = default;

    [[using gnu: always_inline]]
    explicit MatrixX(size_type nrows, size_type ncols)
        : m_data(nrows * ncols), m_nrows{nrows}, m_ncols{ncols} {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    explicit MatrixX(size_type nrows, size_type ncols, const_reference value)
        : m_data(nrows * ncols, value), m_nrows{nrows}, m_ncols{ncols} {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    explicit MatrixX(MatrixView<value_type, kOrder> matrix_view)
        : m_data(matrix_view.nrows() * matrix_view.ncols()), m_nrows{matrix_view.nrows()},
          m_ncols{matrix_view.ncols()} {
        m_data.shrink_to_fit();
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                operator[](i, j) = matrix_view[i, j];
            }
        }
    }

    [[using gnu: always_inline]]
    explicit MatrixX(MatrixView<const value_type, kOrder> matrix_view)
        : m_data(matrix_view.nrows() * matrix_view.ncols()), m_nrows{matrix_view.nrows()},
          m_ncols{matrix_view.ncols()} {
        m_data.shrink_to_fit();
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                operator[](i, j) = matrix_view[i, j];
            }
        }
    }

    [[using gnu: pure, always_inline]] operator MatrixView<value_type, kOrder>() noexcept {
        return MatrixView<value_type, kOrder>(m_nrows, m_ncols, m_data.data());
    }

    [[using gnu: pure,
      always_inline]] operator MatrixView<const value_type, kOrder>() const noexcept {
        return MatrixView<const value_type, kOrder>(m_nrows, m_ncols, m_data.data());
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    auto nrows() const noexcept -> size_type {
        return m_nrows;
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    auto ncols() const noexcept -> size_type {
        return m_ncols;
    }

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    static constexpr auto order() noexcept -> MatrixOrder {
        return kOrder;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto size() const noexcept -> size_type {
        return m_nrows * m_ncols;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto stride() const noexcept -> size_type {
        return kOrder == MatrixOrder::COL_MAJOR ? m_nrows : m_ncols;
    }

    [[using gnu: pure, always_inline]]
    auto data() noexcept -> pointer {
        return m_data.data();
    }
    [[using gnu: pure, always_inline]]
    auto data() const noexcept -> const_pointer {
        return m_data.data();
    }

    [[using gnu: always_inline]]
    auto reshape(size_type nrows, size_type ncols) -> void {
        m_data.resize(nrows * ncols);
        m_data.shrink_to_fit();
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto resize(size_type size) -> void {
        m_data.resize(size);
        m_data.shrink_to_fit();
        if constexpr (kOrder == MatrixOrder::COL_MAJOR) {
            m_nrows = size;
            m_ncols = 1;
        }
        if constexpr (kOrder == MatrixOrder::ROW_MAJOR) {
            m_nrows = 1;
            m_ncols = size;
        }
        return;
    }

    [[using gnu: always_inline]]
    auto assign(size_type nrows, size_type ncols, const_reference value) -> void {
        m_data.assign(nrows * ncols, value);
        m_data.shrink_to_fit();
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto assign(size_type n, const_reference value) -> void {
        m_data.assign(n, value);
        m_data.shrink_to_fit();
        if constexpr (kOrder == MatrixOrder::COL_MAJOR) {
            m_nrows = n;
            m_ncols = 1;
        }
        if constexpr (kOrder == MatrixOrder::ROW_MAJOR) {
            m_nrows = 1;
            m_ncols = n;
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    auto fill(const_reference value) noexcept -> void {
        std::fill(data(), data() + size(), value);
        return;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto empty() noexcept -> bool {
        return size() == 0;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator[](size_type row, size_type col) noexcept -> reference {
        return m_data[offset(row, col)];
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator[](size_type row, size_type col) const noexcept -> const_reference {
        return m_data[offset(row, col)];
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator[](size_type i) noexcept -> reference {
        return m_data[offset(i)];
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator[](size_type i) const noexcept -> const_reference {
        return m_data()[offset(i)];
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto coeff(size_type row, size_type col) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (row >= m_nrows || col >= m_ncols) {
            throw std::out_of_range("Matrix index out of range");
        }
#endif
        return operator[](row, col);
    }

    [[using gnu: always_inline, hot]]
    auto updateCoeff(size_type row, size_type col, const_reference value) noexcept(
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
    auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Matrix index out of range.");
        }
#endif
        return operator[](i);
    }

    [[using gnu: always_inline, hot]]
    auto updateCoeff(size_type i, const_reference value) const noexcept(!BOYLE_CHECK_PARAMS)
        -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Matrix index out of range.");
        }
#endif
        operator[](i) = value;
        return;
    }

    [[using gnu: pure, always_inline]]
    auto view(
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
    auto view(std::pair<size_type, size_type> row_range, std::pair<size_type, size_type> col_range)
        const noexcept -> MatrixView<const value_type, kOrder> {
        row_range.second = std::min(row_range.second, m_nrows);
        col_range.second = std::min(col_range.second, m_ncols);
        return MatrixView<const value_type, kOrder>(
            data() + offset(row_range.first, col_range.first), row_range.second - row_range.first,
            col_range.second - col_range.first, stride()

        );
    }

    [[using gnu: pure, always_inline]]
    auto view(std::pair<size_type, size_type> range) noexcept(!BOYLE_CHECK_PARAMS)
        -> MatrixView<value_type, kOrder> {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != 1 && m_ncols != 1) {
            throw std::invalid_argument("Vector slice only applies on a vector.");
        }
#endif
        return view(
            m_ncols == 1 ? range : std::pair<size_type, size_type>{0, 1},
            m_ncols == 1 ? std::pair<size_type, size_type>{0, 1} : range
        );
    }

    [[using gnu: pure, always_inline]]
    auto view(std::pair<size_type, size_type> range) const noexcept(!BOYLE_CHECK_PARAMS)
        -> MatrixView<const value_type, kOrder> {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != 1 && m_ncols != 1) {
            throw std::invalid_argument("Vector slice only applies on a vector.");
        }
#endif
        return view(
            m_ncols == 1 ? range : std::pair<size_type, size_type>{0, 1},
            m_ncols == 1 ? std::pair<size_type, size_type>{0, 1} : range
        );
    }

    [[using gnu: pure, always_inline]]
    auto view() noexcept -> MatrixView<value_type, kOrder> {
        return view(
            {0, std::numeric_limits<value_type>::max()}, {0, std::numeric_limits<value_type>::max()}
        );
    }

    [[using gnu: pure, always_inline]]
    auto view() const noexcept -> MatrixView<const value_type, kOrder> {
        return view(
            {0, std::numeric_limits<value_type>::max()}, {0, std::numeric_limits<value_type>::max()}
        );
    }

    [[using gnu: always_inline]]
    auto setIdentity() noexcept -> void {
        const size_type n{std::min(m_nrows, m_ncols)};
        fill(0.0);
        for (size_type i{0}; i < n; ++i) {
            operator[](i, i) = 1.0;
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    auto operator+=(const MatrixX& obj) noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
#if BOYLE_USE_BLAS_LAPACK
        constexpr value_type alpha(1.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_saxpy(size(), alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, double>) {
            cblas_daxpy(size(), alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_caxpy(size(), &alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zaxpy(size(), &alpha, obj.data(), 1, data(), 1);
        }
#else
        std::transform(data(), data() + size(), obj.data(), data(), std::plus<value_type>());
#endif
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator+(const MatrixX& obj) const& noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        MatrixX result{*this};
        result += obj;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator+(const MatrixX& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator-=(const MatrixX& obj) noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
#if BOYLE_USE_BLAS_LAPACK
        constexpr value_type alpha(-1.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_saxpy(size(), alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, double>) {
            cblas_daxpy(size(), alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_caxpy(size(), &alpha, obj.data(), 1, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zaxpy(size(), &alpha, obj.data(), 1, data(), 1);
        }
#else
        std::transform(data(), data() + size(), obj.data(), data(), std::minus<value_type>());
#endif
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator-(const MatrixX& obj) const& noexcept -> MatrixX {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        MatrixX result{*this};
        result -= obj;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator-(const MatrixX& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator*=(const ScalarArithmetic auto& fac) noexcept -> MatrixX& {
        const value_type alpha(fac);
#if BOYLE_USE_BLAS_LAPACK
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sscal(size(), alpha, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, double>) {
            cblas_dscal(size(), alpha, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cscal(size(), &alpha, data(), 1);
        }
        if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zscal(size(), &alpha, data(), 1);
        }
#else
        std::transform(
            data(), data() + size(), data(),
            [alpha](const value_type& x) noexcept -> value_type { return x * alpha; }
        );
#endif
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator*(const ScalarArithmetic auto& fac) const& noexcept -> MatrixX {
        MatrixX result{*this};
        result *= fac;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator*(const ScalarArithmetic auto& fac) && noexcept -> MatrixX&& {
        operator*=(fac);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator/=(const ScalarArithmetic auto& den) noexcept -> MatrixX& {
        const auto fac{decltype(den)(1.0) / den};
        operator*=(fac);
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator/(ScalarArithmetic auto den) const& noexcept -> MatrixX {
        MatrixX result{*this};
        result /= den;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator/(ScalarArithmetic auto den) && noexcept -> MatrixX&& {
        operator/=(den);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator-() const& noexcept -> MatrixX {
        Matrix result{*this};
        result *= -1.0;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator-() && noexcept -> MatrixX&& {
        operator*=(-1.0);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, hot]]
    auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type> {
        detail::DenseNormTraitT<value_type> result(0.0);
        if (m_nrows == 1 && m_ncols == 1) {
            result = std::abs(m_data[0]);
        }
#ifdef BOYLE_USE_BLAS_LAPACK
        if (m_nrows == 1 || m_ncols == 1) {
            if constexpr (std::is_same_v<value_type, float>) {
                result = cblas_snrm2(size(), data(), 1);
            }
            if constexpr (std::is_same_v<value_type, double>) {
                result = cblas_dnrm2(size(), data(), 1);
            }
            if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                result = cblas_scnrm2(size(), data(), 1);
            }
            if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                result = cblas_dznrm2(size(), data(), 1);
            }
        }
#else
        if (m_nrows == 1 || m_ncols == 1) {
            result = std::transform_reduce(
                data(), data() + size(), 0.0, std::plus<detail::DenseNormTraitT<value_type>>{},
                [](const value_type& a) noexcept -> detail::DenseNormTraitT<value_type> {
                    return std::norm(a);
                }
            );
            result = std::sqrt(result);
        }
#endif
        return result;
    }

    [[using gnu: pure, always_inline, hot]]
    auto euclideanSqr() const noexcept -> detail::DenseNormTraitT<value_type> {
        detail::DenseNormTraitT<value_type> result(0.0);
        if (m_nrows == 1 && m_ncols == 1) {
            result = std::norm(m_data[0]);
        }
#ifdef BOYLE_USE_BLAS_LAPACK
        if (m_nrows == 1 || m_ncols == 1) {
            if constexpr (std::is_same_v<value_type, float>) {
                result = cblas_sdot(size(), data(), 1, data(), 1);
            }
            if constexpr (std::is_same_v<value_type, double>) {
                result = cblas_ddot(size(), data(), 1, data(), 1);
            }
            if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                result = cblas_cdotc(size(), data(), 1, data(), 1);
            }
            if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                result = cblas_zdotc(size(), data(), 1, data(), 1);
            }
        }
#else
        if (m_nrows == 1 || m_ncols == 1) {
            result = std::transform_reduce(
                data(), data() + size(), 0.0, std::plus<detail::DenseNormTraitT<value_type>>{},
                [](const value_type& a) constexpr noexcept -> detail::DenseNormTraitT<value_type> {
                    return std::norm(a);
                }
            );
        }
#endif
        return result;
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
                test, std::abs(operator[](i) - obj[i]) / std::max(std::abs(operator[](i), 1.0))
            );
        }
        return test < tol;
    }

    [[using gnu: pure, always_inline]]
    auto transposed() const noexcept -> MatrixX {
        MatrixX result(m_ncols, m_nrows);
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                result[j, i] = m_data[offset(i, j)];
            }
        }
        return result;
    }

    [[using gnu: always_inline]]
    auto selfConjugated() noexcept -> MatrixX& {
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            std::transform(
                data(), data() + size(), data(),
                [](const value_type& a) noexcept -> value_type { return std::conj(a); }
            );
        }
        return *this;
    }

    [[using gnu: pure, always_inline]]
    auto conjugated() const& noexcept -> MatrixX {
        MatrixX result{*this};
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            result.selfConjugated();
        }
        return result;
    }

    [[using gnu: always_inline]]
    auto conjugated() && noexcept -> MatrixX&& {
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            selfConjugated();
        }
        return std::move(*this);
    }

    [[using gnu: pure, always_inline]]
    auto adjoint() const noexcept -> MatrixX {
        return transposed().selfConjugated();
    }

    [[using gnu: pure, always_inline, hot]]
    auto dot(const MatrixX& obj) const noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ncols != obj.nrows()) {
            throw std::invalid_argument("MatrixX::dot: incompatible dimensions");
        }
#endif
        const std::size_t outer_size = obj.ncols();
        MatrixX result(m_nrows, outer_size, 0.0);
#ifdef BOYLE_USE_BLAS_LAPACK
        constexpr CBLAS_ORDER order =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? CblasColMajor : CblasRowMajor);
        const blasint lda = (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : m_ncols);
        const blasint ldb =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_ncols : outer_size);
        const blasint ldc =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : outer_size);
        constexpr value_type alpha(1.0), beta(0.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, alpha, data(), lda,
                obj.data(), ldb, beta, result.data(), ldc
            );
        }
        if constexpr (std::is_same_v<value_type, double>) {
            cblas_dgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, alpha, data(), lda,
                obj.data(), ldb, beta, result.data(), ldc
            );
        }
        if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, &alpha, data(),
                lda, obj.data(), ldb, &beta, result.data(), ldc
            );
        }
        if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, &alpha, data(),
                lda, obj.data(), ldb, &beta, result.data(), ldc
            );
        }
#else
        for (std::size_t i{0}; i < m_nrows; ++i) {
            for (std::size_t j{0}; j < outer_size; ++j) {
                for (std::size_t k{0}; k < m_ncols; ++k) {
                    result[i, j] += operator[](i, k) * obj[k, j];
                }
            }
        }
#endif
        return result;
    }

    [[using gnu: pure, always_inline, hot]]
    auto dot(const VectorX<value_type, allocator_type>& obj) const noexcept(!BOYLE_CHECK_PARAMS)
        -> VectorX<value_type, allocator_type> {
        VectorX<value_type, allocator_type> result(m_nrows, 0.0);
#if BOYLE_CHECK_PARAMS == 1
        if (m_ncols != obj.size()) {
            throw std::invalid_argument("MatrixX::dot: incompatible dimensions");
        }
#endif
#ifdef BOYLE_USE_BLAS_LAPACK
        constexpr CBLAS_ORDER order =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? CblasColMajor : CblasRowMajor);
        blasint lda = (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : m_ncols);
        constexpr value_type alpha(1.0), beta(0.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sgemv(
                order, CblasNoTrans, m_nrows, m_ncols, alpha, data(), lda, obj.data(), 1, beta,
                result.data(), 1
            );
        }
        if constexpr (std::is_same_v<value_type, double>) {
            cblas_dgemv(
                order, CblasNoTrans, m_nrows, m_ncols, alpha, data(), lda, obj.data(), 1, beta,
                result.data(), 1
            );
        }
        if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cgemv(
                order, CblasNoTrans, m_nrows, m_ncols, &alpha, data(), lda, obj.data(), 1, &beta,
                result.data(), 1
            );
        }
        if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zgemv(
                order, CblasNoTrans, m_nrows, m_ncols, &alpha, data(), lda, obj.data(), 1, &beta,
                result.data(), 1
            );
        }
#else
        for (std::size_t i{0}; i < m_nrows; ++i) {
            for (std::size_t j{0}; j < m_ncols; ++j) {
                result[i] += operator[](i, j) * obj[j];
            }
        }
#endif
        return result;
    }

  private:
    [[using gnu: pure, always_inline, leaf, hot]] [[nodiscard]]
    auto offset(size_type row, size_type col) const noexcept -> size_type {
        if constexpr (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR) {
            return (col * m_nrows) + row;
        }
        if constexpr (kOrder == ::boyle::math::MatrixOrder::ROW_MAJOR) {
            return (row * m_ncols) + col;
        }
    }

    [[using gnu: pure, always_inline, leaf, hot]] [[nodiscard]]
    auto offset(size_type i) const noexcept -> size_type {
        if constexpr (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR) {
            return i / m_nrows * m_nrows + i % m_nrows;
        }
        if constexpr (kOrder == ::boyle::math::MatrixOrder::ROW_MAJOR) {
            return i / m_ncols * m_ncols + i % m_ncols;
        }
    }

    std::vector<value_type, allocator_type> m_data{};
    size_type m_nrows{0}, m_ncols{0};
};

template <ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
[[using gnu: pure, always_inline, hot]]
inline auto operator*(
    const ScalarArithmetic auto& fac, const MatrixX<Scalar, Order, Alloc>& obj
) noexcept -> MatrixX<Scalar, Order, Alloc> {
    return obj * fac;
}

template <ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
[[using gnu: always_inline, hot]]
inline auto operator*(
    const ScalarArithmetic auto& fac, MatrixX<Scalar, Order, Alloc>&& obj
) noexcept -> MatrixX<Scalar, Order, Alloc>&& {
    return std::move(obj) * fac;
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

template <boyle::math::ScalarArithmetic Scalar, boyle::math::MatrixOrder Order, typename Alloc>
[[using gnu: pure, always_inline]]
inline auto abs(const boyle::math::MatrixX<Scalar, Order, Alloc>& matrix) noexcept
    -> boyle::math::detail::DenseNormTrait<boyle::math::MatrixX<Scalar, Order, Alloc>> {
    boyle::math::detail::DenseNormTraitT<boyle::math::MatrixX<Scalar, Order, Alloc>> result(matrix);
    if constexpr (boyle::math::isComplexArithmeticV<Scalar>) {
        std::transform(
            result.data(), result.data() + result.size(), result.data(),
            [](const Scalar& x) noexcept -> boyle::math::detail::DenseNormTraitT<Scalar> {
                return std::abs(x);
            }
        );
    }
    return result;
}

template <boyle::math::ScalarArithmetic Scalar, boyle::math::MatrixOrder Order, typename Alloc>
[[using gnu: pure, always_inline]]
inline auto conj(const boyle::math::MatrixX<Scalar, Order, Alloc>& matrix) noexcept
    -> boyle::math::MatrixX<Scalar, Order, Alloc> {
    return matrix.conjugated();
}

template <boyle::math::ScalarArithmetic Scalar, boyle::math::MatrixOrder Order, typename Alloc>
[[using gnu: always_inline]]
inline auto conj(boyle::math::MatrixX<Scalar, Order, Alloc>&& matrix) noexcept
    -> boyle::math::MatrixX<Scalar, Order, Alloc>&& {
    return std::move(matrix).conjugated();
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
