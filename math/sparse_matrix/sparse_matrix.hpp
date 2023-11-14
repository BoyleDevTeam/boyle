/**
 * @file sparse_matrix.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-19
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <unordered_map>
#include <vector>

#include "common/utils/macros.hpp"
#include "math/type_traits.hpp"

namespace tiny_pnc {
namespace math {

template <typename DerivedMatrix, typename Scalar = double, typename Index = int>
class [[nodiscard]] SparseMatrix {
    static_assert(
        std::is_arithmetic_v<Scalar> || isVecArithmeticV<Scalar>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_integral_v<Index>, "The loaded type must be a floating-point type.");

  public:
    ENABLE_IMPLICIT_CONSTRUCTORS(SparseMatrix);
    virtual ~SparseMatrix() noexcept = default;

    [[using gnu: pure, always_inline]]
    std::size_t nrows() const noexcept {
        return underlying()->nrows();
    }

    [[using gnu: pure, always_inline]]
    std::size_t ncols() const noexcept {
        return underlying()->ncols();
    }

    [[using gnu: pure, always_inline]]
    std::size_t nnzs() const noexcept {
        return underlying()->nnzs();
    }

    [[using gnu: always_inline]]
    void resize(std::size_t nrows, std::size_t ncols) noexcept {
        underlying()->resize(nrows, ncols);
        return;
    }

    [[using gnu: pure, always_inline]]
    Scalar coeff(Index row, Index col) const noexcept {
        return underlying()->coeff(row, col);
    }

    [[using gnu: pure, always_inline]]
    Scalar
    operator()(Index row, Index col) const noexcept {
        return coeff(row, col);
    }

  private:
    [[using gnu: always_inline]]
    DerivedMatrix* underlying() noexcept {
        return static_cast<DerivedMatrix*>(this);
    }

    [[using gnu: always_inline]]
    const DerivedMatrix* underlying() const noexcept {
        return static_cast<const DerivedMatrix*>(this);
    }
};

} // namespace math
} // namespace tiny_pnc
