/**
 * @file csr_matrix.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-20
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <map>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "spdlog/spdlog.h"

#include "common/utils/macros.hpp"
#include "math/sparse_matrix/sparse_matrix.hpp"

namespace tiny_pnc {
namespace math {

template <typename T, typename Index>
class CooMatrix;

template <typename T, typename Index>
class LilMatrix;

template <typename T = double, typename Index = int>
class [[nodiscard]] CsrMatrix final : public SparseMatrix<CsrMatrix<T, Index>, T, Index> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_integral_v<Index>, "The loaded type must be a floating-point type.");

    friend class boost::serialization::access;

  public:
    [[using gnu: always_inline]] explicit CsrMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : ncols_{ncols}, values_{}, inner_indices_{}, outer_indices_(nrows + 1, 0) {}

    ENABLE_IMPLICIT_CONSTRUCTORS(CsrMatrix);

    ~CsrMatrix() noexcept override = default;

    [[using gnu: flatten, leaf]] CsrMatrix(const CooMatrix<T, Index>& coo_matrix) noexcept
        : ncols_{coo_matrix.ncols_}, values_{}, inner_indices_{}, outer_indices_{} {
        std::map<
            typename CooMatrix<T, Index>::IndexPair, std::size_t,
            typename CooMatrix<T, Index>::IndexPairRowMajorCompare>
            offsets_map{coo_matrix.offsets_map_.cbegin(), coo_matrix.offsets_map_.cend()};
        const std::size_t nnzs = coo_matrix.nnzs();
        values_.reserve(nnzs);
        inner_indices_.reserve(nnzs);
        outer_indices_.reserve(coo_matrix.nrows_ + 1);
        outer_indices_.push_back(0);
        for (const auto& [index_pair, offset] : offsets_map) {
            values_.push_back(coo_matrix.values_[offset]);
            inner_indices_.push_back(index_pair.col);
            const std::size_t diff_row = index_pair.row - outer_indices_.size() + 1;
            for (std::size_t i{0}; i < diff_row; ++i) {
                outer_indices_.push_back(values_.size() - 1);
            }
        }
        for (std::size_t i = outer_indices_.size(); i < coo_matrix.nrows_ + 1; ++i) {
            outer_indices_.push_back(nnzs);
        }
    }

    [[using gnu: always_inline]] CsrMatrix(const LilMatrix<T, Index>& lil_matrix) noexcept {
        CooMatrix<T, Index> coo_matrix{lil_matrix};
        *this = coo_matrix;
    }

    [[using gnu: pure, always_inline]]
    std::size_t nrows() const noexcept {
        return outer_indices_.size() - 1;
    }

    [[using gnu: pure, always_inline]]
    std::size_t ncols() const noexcept {
        return ncols_;
    }

    [[using gnu: pure, always_inline]]
    std::size_t nnzs() const noexcept {
        return values_.size();
    }

    [[using gnu: flatten, leaf]]
    void resize(std::size_t nrows, std::size_t ncols) noexcept {
        if (nrows >= outer_indices_.size() - 1) {
            for (std::size_t i = outer_indices_.size() - 1; i < nrows; ++i) {
                outer_indices_.push_back(outer_indices_.back());
            }
        } else {
            outer_indices_.resize(nrows + 1);
            values_.resize(outer_indices_.back());
            inner_indices_.resize(outer_indices_.back());
        }
        const std::size_t values_size{values_.size()};
        if (ncols < ncols_) {
            std::size_t j{0}, k{0};
            ncols_ = ncols;
            for (std::size_t i{0}; i < values_size; ++i) {
                if (inner_indices_[i] >= static_cast<Index>(ncols_)) {
                    continue;
                }
                if (j != i) {
                    values_[j] = values_[i];
                    inner_indices_[j] = inner_indices_[i];
                    if (i == static_cast<std::size_t>(outer_indices_[k])) {
                        outer_indices_[k] = j;
                    }
                }
                ++j;
                if (i == static_cast<std::size_t>(outer_indices_[k])) {
                    ++k;
                }
            }
            outer_indices_[ncols] = j;
            values_.resize(j);
            inner_indices_.resize(j);
            outer_indices_.shrink_to_fit();
            values_.shrink_to_fit();
            inner_indices_.shrink_to_fit();
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    T coeff(Index row, Index col) const noexcept {
        if (row >= static_cast<Index>(outer_indices_.size() - 1) ||
            col >= static_cast<Index>(ncols_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, outer_indices_.size() - 1, ncols_
            );
            return 0.0;
        }
        if (const auto search = std::find(
                inner_indices_.cbegin() + outer_indices_[row],
                inner_indices_.cbegin() + outer_indices_[row + 1], col
            );
            search != inner_indices_.cbegin() + outer_indices_[row + 1]) {
            return values_[search - inner_indices_.cbegin()];
        }
        return 0.0;
    }

    [[using gnu: hot]]
    void updateCoeff(Index row, Index col, T value) noexcept {
        if (row >= static_cast<Index>(outer_indices_.size() - 1) ||
            col >= static_cast<Index>(ncols_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, outer_indices_.size() - 1, ncols_
            );
            return;
        }
        if (const auto search = std::find(
                inner_indices_.cbegin() + outer_indices_[row],
                inner_indices_.cbegin() + outer_indices_[row + 1], col
            );
            search != inner_indices_.cbegin() + outer_indices_[row + 1]) {
            const std::size_t inner_offset = search - inner_indices_.cbegin();
            if (value != 0.0) {
                values_[inner_offset] = value;
            } else {
                values_.erase(values_.begin() + inner_offset);
                inner_indices_.erase(inner_indices_.begin() + inner_offset);
                for (auto outer_it = outer_indices_.begin() + row + 1;
                     outer_it != outer_indices_.end(); ++outer_it) {
                    --(*outer_it);
                }
            }
            return;
        } else {
            if (value != 0.0) {
                const std::size_t inner_offset =
                    std::upper_bound(
                        inner_indices_.cbegin() + outer_indices_[row],
                        inner_indices_.cbegin() + outer_indices_[row + 1], col
                    ) -
                    inner_indices_.cbegin();
                values_.insert(values_.begin() + inner_offset, value);
                inner_indices_.insert(inner_indices_.begin() + inner_offset, col);
                for (auto outer_it = outer_indices_.begin() + row + 1;
                     outer_it != outer_indices_.end(); ++outer_it) {
                    ++(*outer_it);
                }
            }
        }
        return;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& values() const noexcept {
        return values_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<Index>& innerIndices() const noexcept {
        return inner_indices_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<Index>& outerIndices() const noexcept {
        return outer_indices_;
    }

    [[using gnu: always_inline]]
    void clear() noexcept {
        values_.clear();
        inner_indices_.clear();
        std::fill(outer_indices_.begin(), outer_indices_.end(), 0);
        return;
    }

  private:
    template <typename Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) {
        ar& ncols_;
        ar& values_;
        ar& inner_indices_;
        ar& outer_indices_;
        return;
    }

    std::size_t ncols_;
    std::vector<T> values_;
    std::vector<Index> inner_indices_;
    std::vector<Index> outer_indices_;
};

} // namespace math
} // namespace tiny_pnc
