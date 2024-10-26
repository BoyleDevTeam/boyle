/**
 * @file partial_pivot.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-01-22
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#ifndef __MATH_PARTIAL_PIVOT_HPP__
#define __MATH_PARTIAL_PIVOT_HPP__

#include <cmath>
#include <numeric>

#include "poissonsoft/math/dense/detail/dense_norm_result.hpp"
#include "poissonsoft/math/dense/detail/dense_partial_pivot_result.hpp"

namespace poissonsoft::math {

template <typename T>
inline auto partialPivot(const T& matrix) noexcept(!POISSONSOFT_CHECK_PARAMS)
    -> detail::DensePartialPivotTraitT<T> {
    using matrix_type = T;
    using value_type = typename matrix_type::value_type;
    using size_type = typename matrix_type::size_type;

    const size_type n{matrix.nrows()};
    detail::DensePartialPivotTraitT<matrix_type> ipiv(n);
    std::iota(ipiv.begin(), ipiv.end(), 0);
    for (size_type k{0}; k < n; ++k) {
        detail::DenseNormTraitT<value_type> max{std::abs(matrix(ipiv[k], k))};
        size_type imax{k};
        for (size_type i{k + 1}; i < n; ++i) {
            detail::DenseNormTraitT<value_type> temp{std::abs(matrix(ipiv[i], k))};
            if (temp > max) {
                max = temp;
                imax = i;
            }
        }
        if (imax != k) {
            std::swap(ipiv[k], ipiv[imax]);
        }
    }
    return ipiv;
}

} // namespace poissonsoft::math

#endif
