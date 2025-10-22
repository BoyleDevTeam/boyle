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
#include <iomanip>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/complex.hpp"
#include "boost/serialization/vector.hpp"

#ifdef BOYLE_USE_BLAS_LAPACK
#include "cblas.h"
#include "lapacke.h"
#endif

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/dense/matrix_view.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
class MatrixX final {
    friend class boost::serialization::access;

  public:
    using value_type = typename DenseTraits<MatrixX>::value_type;
    using reference = typename DenseTraits<MatrixX>::reference;
    using const_reference = typename DenseTraits<MatrixX>::const_reference;
    using pointer = typename DenseTraits<MatrixX>::pointer;
    using const_pointer = typename DenseTraits<MatrixX>::const_pointer;
    using size_type = typename DenseTraits<MatrixX>::size_type;
    using difference_type = typename DenseTraits<MatrixX>::difference_type;
    using allocator_type = typename DenseTraits<MatrixX>::allocator_type;

    static constexpr MatrixOrder kOrder = Order;

    [[using gnu: always_inline]]
    explicit MatrixX(const allocator_type& alloc = {}) noexcept
        : m_data(alloc) {}
    MatrixX(const MatrixX& other) = default;
    auto operator=(const MatrixX& other) -> MatrixX& = default;
    MatrixX(MatrixX&& other) noexcept = default;
    auto operator=(MatrixX&& other) noexcept -> MatrixX& = default;
    ~MatrixX() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_data.get_allocator();
    }

    [[using gnu: always_inline]]
    explicit MatrixX(size_type nrows, size_type ncols, const allocator_type& alloc = {})
        : m_data(nrows * ncols, alloc), m_nrows{nrows}, m_ncols{ncols} {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    explicit MatrixX(
        size_type nrows, size_type ncols, const_reference value, const allocator_type& alloc = {}
    )
        : m_data(nrows * ncols, value, alloc), m_nrows{nrows}, m_ncols{ncols} {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    explicit MatrixX(MatrixView<value_type, kOrder> matrix_view, const allocator_type& alloc = {})
        : m_data(matrix_view.nrows() * matrix_view.ncols(), alloc), m_nrows{matrix_view.nrows()},
          m_ncols{matrix_view.ncols()} {
        m_data.shrink_to_fit();
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                operator[](i, j) = matrix_view[i, j];
            }
        }
    }

    [[using gnu: always_inline]]
    explicit MatrixX(
        MatrixView<const value_type, kOrder> matrix_view, const allocator_type& alloc = {}
    )
        : m_data(matrix_view.nrows() * matrix_view.ncols(), alloc), m_nrows{matrix_view.nrows()},
          m_ncols{matrix_view.ncols()} {
        m_data.shrink_to_fit();
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
                operator[](i, j) = matrix_view[i, j];
            }
        }
    }

    [[using gnu: pure, always_inline]] operator MatrixView<value_type, kOrder>() noexcept {
        return MatrixView<value_type, kOrder>(m_data.data(), m_nrows, m_ncols);
    }

    [[using gnu: pure, always_inline]] operator MatrixView<
        const value_type, kOrder>() const noexcept {
        return MatrixView<const value_type, kOrder>(m_data.data(), m_nrows, m_ncols);
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
    auto resize(size_type nrows, size_type ncols) -> void {
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
        return m_data[offset(i)];
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto coeff(size_type row, size_type col) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
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
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            throw std::out_of_range("Matrix index out of range");
        }
#endif
        operator[](row, col) = value;
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) [[unlikely]] {
            throw std::out_of_range("Matrix index out of range.");
        }
#endif
        return operator[](i);
    }

    [[using gnu: always_inline, hot]]
    auto updateCoeff(size_type i, const_reference value) const noexcept(!BOYLE_CHECK_PARAMS)
        -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) [[unlikely]] {
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
    auto view(
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
    auto view(std::pair<size_type, size_type> range) noexcept(!BOYLE_CHECK_PARAMS)
        -> MatrixView<value_type, kOrder> {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != 1 && m_ncols != 1) [[unlikely]] {
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
        if (m_nrows != 1 && m_ncols != 1) [[unlikely]] {
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
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
#ifdef BOYLE_USE_BLAS_LAPACK
        [[maybe_unused]] constexpr value_type alpha(1.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_saxpy(size(), alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, double>) {
            cblas_daxpy(size(), alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_caxpy(size(), &alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zaxpy(size(), &alpha, obj.data(), 1, data(), 1);
        } else {
            std::transform(data(), data() + size(), obj.data(), data(), std::plus<value_type>());
        }
#else
        std::transform(data(), data() + size(), obj.data(), data(), std::plus<value_type>());
#endif
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator+(const MatrixX& obj) const& noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        MatrixX result{*this};
        result += obj;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator+(MatrixX&& obj) const& noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        obj += *this;
        return std::move(obj);
    }

    [[using gnu: always_inline, hot]]
    auto operator+(const MatrixX& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator+(MatrixX&& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator-=(const MatrixX& obj) noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
#ifdef BOYLE_USE_BLAS_LAPACK
        [[maybe_unused]] constexpr value_type alpha(-1.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_saxpy(size(), alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, double>) {
            cblas_daxpy(size(), alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_caxpy(size(), &alpha, obj.data(), 1, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zaxpy(size(), &alpha, obj.data(), 1, data(), 1);
        } else {
            std::transform(data(), data() + size(), obj.data(), data(), std::minus<value_type>());
        }
#else
        std::transform(data(), data() + size(), obj.data(), data(), std::minus<value_type>());
#endif
        return *this;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator-(const MatrixX& obj) const& noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        MatrixX result{*this};
        result -= obj;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator-(MatrixX&& obj) const& noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        obj *= -1.0;
        obj += *this;
        return std::move(obj);
    }

    [[using gnu: always_inline, hot]]
    auto operator-(const MatrixX& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator-(MatrixX&& obj) && noexcept(!BOYLE_CHECK_PARAMS) -> MatrixX&& {
#if BOYLE_CHECK_PARAMS == 1
        if (m_nrows != obj.m_nrows || m_ncols != obj.m_ncols) [[unlikely]] {
            throw std::invalid_argument("Matrix dimensions do not match.");
        }
#endif
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    auto operator*=(const ScalarArithmetic auto& fac) noexcept -> MatrixX& {
        const auto alpha{static_cast<value_type>(fac)};
#ifdef BOYLE_USE_BLAS_LAPACK
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sscal(size(), alpha, data(), 1);
        } else if constexpr (std::is_same_v<value_type, double>) {
            cblas_dscal(size(), alpha, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cscal(size(), &alpha, data(), 1);
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zscal(size(), &alpha, data(), 1);
        } else {
            std::transform(
                data(), data() + size(), data(),
                [alpha](const value_type& x) constexpr noexcept -> value_type { return x * alpha; }
            );
        }
#else
        std::transform(
            data(), data() + size(), data(),
            [alpha](const value_type& x) constexpr noexcept -> value_type { return x * alpha; }
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
        if constexpr (std::is_integral_v<value_type>) {
            const auto alpha{static_cast<value_type>(den)};
            std::transform(
                data(), data() + size(), data(),
                [alpha](const value_type& x) constexpr noexcept -> value_type { return x / alpha; }
            );
        } else {
            const auto alpha{static_cast<value_type>(1.0 / den)};
            operator*=(alpha);
        }
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
        MatrixX result{*this};
        result *= -1.0;
        return result;
    }

    [[using gnu: always_inline, hot]]
    auto operator-() && noexcept -> MatrixX&& {
        operator*=(-1.0);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator==(const MatrixX& other) const noexcept -> bool {
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
    auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type> {
        detail::DenseNormTraitT<value_type> result(0.0);
        if (m_nrows == 1 && m_ncols == 1) {
            result = std::abs(m_data[0]);
        }
#ifdef BOYLE_USE_BLAS_LAPACK
        if (m_nrows == 1 || m_ncols == 1) {
            if constexpr (std::is_same_v<value_type, float>) {
                result = cblas_snrm2(size(), data(), 1);
            } else if constexpr (std::is_same_v<value_type, double>) {
                result = cblas_dnrm2(size(), data(), 1);
            } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                result = cblas_scnrm2(size(), data(), 1);
            } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                result = cblas_dznrm2(size(), data(), 1);
            } else {
                result = std::transform_reduce(
                    data(), data() + size(), 0.0, std::plus<detail::DenseNormTraitT<value_type>>{},
                    [](const value_type& a) constexpr noexcept
                        -> detail::DenseNormTraitT<value_type> { return std::norm(a); }
                );
                result = std::sqrt(result);
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
            } else if constexpr (std::is_same_v<value_type, double>) {
                result = cblas_ddot(size(), data(), 1, data(), 1);
            } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                result = cblas_cdotc(size(), data(), 1, data(), 1);
            } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                result = cblas_zdotc(size(), data(), 1, data(), 1);
            } else {
                result = std::transform_reduce(
                    data(), data() + size(), 0.0, std::plus<detail::DenseNormTraitT<value_type>>{},
                    [](const value_type& a) constexpr noexcept
                        -> detail::DenseNormTraitT<value_type> { return std::norm(a); }
                );
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
                test, std::abs(operator[](i) - obj[i]) /
                          std::max(std::abs(operator[](i)), static_cast<decltype(test)>(1.0))
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
#ifdef BOYLE_USE_BLAS_LAPACK
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                LAPACKE_clacgv_work(size(), reinterpret_cast<lapack_complex_float*>(data()), 1);
            } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                LAPACKE_zlacgv_work(size(), reinterpret_cast<lapack_complex_double*>(data()), 1);
            } else {
                std::transform(
                    data(), data() + size(), data(),
                    [](const value_type& a) constexpr noexcept -> value_type {
                        return std::conj(a);
                    }
                );
            }
        }
#else
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            std::transform(
                data(), data() + size(), data(),
                [](const value_type& a) constexpr noexcept -> value_type { return std::conj(a); }
            );
        }
#endif
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
        if (m_ncols != obj.nrows()) [[unlikely]] {
            throw std::invalid_argument("MatrixX::dot: incompatible dimensions");
        }
#endif
        const size_type outer_size = obj.ncols();
        MatrixX result(m_nrows, outer_size, 0.0);
#ifdef BOYLE_USE_BLAS_LAPACK
        [[maybe_unused]] constexpr CBLAS_ORDER order =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? CblasColMajor : CblasRowMajor);
        [[maybe_unused]] const blasint lda =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : m_ncols);
        [[maybe_unused]] const blasint ldb =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_ncols : outer_size);
        [[maybe_unused]] const blasint ldc =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : outer_size);
        [[maybe_unused]] constexpr value_type alpha(1.0), beta(0.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, alpha, data(), lda,
                obj.data(), ldb, beta, result.data(), ldc
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            cblas_dgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, alpha, data(), lda,
                obj.data(), ldb, beta, result.data(), ldc
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, &alpha, data(),
                lda, obj.data(), ldb, &beta, result.data(), ldc
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zgemm(
                order, CblasNoTrans, CblasNoTrans, m_nrows, outer_size, m_ncols, &alpha, data(),
                lda, obj.data(), ldb, &beta, result.data(), ldc
            );
        } else {
            for (size_type i{0}; i < m_nrows; ++i) {
                for (size_type j{0}; j < outer_size; ++j) {
                    for (size_type k{0}; k < m_ncols; ++k) {
                        result[i, j] += operator[](i, k) * obj[k, j];
                    }
                }
            }
        }
#else
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < outer_size; ++j) {
                for (size_type k{0}; k < m_ncols; ++k) {
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
        if (m_ncols != obj.size()) [[unlikely]] {
            throw std::invalid_argument("MatrixX::dot: incompatible dimensions");
        }
#endif
#ifdef BOYLE_USE_BLAS_LAPACK
        [[maybe_unused]] constexpr CBLAS_ORDER order =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? CblasColMajor : CblasRowMajor);
        [[maybe_unused]] blasint lda =
            (kOrder == ::boyle::math::MatrixOrder::COL_MAJOR ? m_nrows : m_ncols);
        [[maybe_unused]] constexpr value_type alpha(1.0), beta(0.0);
        if constexpr (std::is_same_v<value_type, float>) {
            cblas_sgemv(
                order, CblasNoTrans, m_nrows, m_ncols, alpha, data(), lda, obj.data(), 1, beta,
                result.data(), 1
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            cblas_dgemv(
                order, CblasNoTrans, m_nrows, m_ncols, alpha, data(), lda, obj.data(), 1, beta,
                result.data(), 1
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            cblas_cgemv(
                order, CblasNoTrans, m_nrows, m_ncols, &alpha, data(), lda, obj.data(), 1, &beta,
                result.data(), 1
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            cblas_zgemv(
                order, CblasNoTrans, m_nrows, m_ncols, &alpha, data(), lda, obj.data(), 1, &beta,
                result.data(), 1
            );
        } else {
            for (size_type i{0}; i < m_nrows; ++i) {
                for (size_type j{0}; j < m_ncols; ++j) {
                    result[i] += operator[](i, j) * obj[j];
                }
            }
        }
#else
        for (size_type i{0}; i < m_nrows; ++i) {
            for (size_type j{0}; j < m_ncols; ++j) {
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

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_data;
        archive & m_nrows;
        archive & m_ncols;
        return;
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
    obj *= fac;
    return std::move(obj);
}

template <typename Char, ScalarArithmetic Scalar, MatrixOrder Order, typename Alloc>
inline auto operator<<(
    std::basic_ostream<Char>& os, const MatrixX<Scalar, Order, Alloc>& matrix
) noexcept -> std::basic_ostream<Char>& {
    using value_type = std::remove_const_t<Scalar>;
    using size_type = typename MatrixX<Scalar, Order, Alloc>::size_type;
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
            [](const Scalar& x) constexpr noexcept -> boyle::math::detail::DenseNormTraitT<Scalar> {
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
    matrix.selfConjugated();
    return std::move(matrix);
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
