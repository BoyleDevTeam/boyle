/**
 * @file sparse_matrix_proxy.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-27
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <utility>

#include "proxy.h"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_nrows, nrows);
PRO_DEF_MEM_DISPATCH(c_ncols, ncols);
PRO_DEF_MEM_DISPATCH(c_nnzs, nnzs);
PRO_DEF_MEM_DISPATCH(m_resize, resize);
PRO_DEF_MEM_DISPATCH(m_reserve, reserve);
PRO_DEF_MEM_DISPATCH(m_clear, clear);
PRO_DEF_MEM_DISPATCH(m_compress, compress);
PRO_DEF_MEM_DISPATCH(c_coeff, coeff);
PRO_DEF_MEM_DISPATCH(m_updateCoeff, updateCoeff);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <::boyle::math::GeneralArithmetic Scalar = double, std::integral Index = int>
class [[nodiscard]] SparseMatrixFacade final : public pro::facade_builder
    ::add_convention<c_nrows, auto() const noexcept -> std::size_t>
    ::add_convention<c_ncols, auto() const noexcept -> std::size_t>
    ::add_convention<c_nnzs, auto() const noexcept -> std::size_t>
    ::add_convention<m_resize, auto(std::size_t, std::size_t) noexcept -> void>
    ::add_convention<m_reserve, auto(std::size_t) noexcept -> void>
    ::add_convention<m_clear, auto() noexcept -> void>
    ::add_convention<m_compress, auto() noexcept -> void>
    ::template add_convention<c_coeff, auto(Index, Index) const noexcept -> Scalar>
    ::template add_convention<m_updateCoeff, auto(Index, Index, Scalar) noexcept -> void>
    ::build {};

// clang-format on

} // namespace detail

template <GeneralArithmetic Scalar = double, std::integral Index = int>
using SparseMatrixProxy = pro::proxy<detail::SparseMatrixFacade<Scalar, Index>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeSparseMatrixProxy(T sparse_matrix
) -> SparseMatrixProxy<typename T::value_type, typename T::index_type> {
    return pro::make_proxy<
        detail::SparseMatrixFacade<typename T::value_type, typename T::index_type>>(
        std::move(sparse_matrix)
    );
}

} // namespace boyle::math
