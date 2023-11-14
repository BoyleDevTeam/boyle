/**
 * @file coo_matrix.hpp
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

#include <unordered_map>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/vector.hpp"
#include "spdlog/spdlog.h"

#include "common/utils/macros.hpp"
#include "math/sparse_matrix/sparse_matrix.hpp"

namespace tiny_pnc {
namespace math {

template <typename T, typename Index>
class LilMatrix;

template <typename T, typename Index>
class CscMatrix;

template <typename T, typename Index>
class CsrMatrix;

template <typename T = double, typename Index = int>
class [[nodiscard]] CooMatrix final : public SparseMatrix<CooMatrix<T, Index>, T, Index> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_integral_v<Index>, "The loaded type must be a floating-point type.");

    friend class boost::serialization::access;
    friend class LilMatrix<T, Index>;
    friend class CscMatrix<T, Index>;
    friend class CsrMatrix<T, Index>;

  public:
    [[using gnu: always_inline]] explicit CooMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : nrows_{nrows}, ncols_{ncols}, offsets_map_{}, values_{}, row_indices_{}, col_indices_{} {}

    ENABLE_IMPLICIT_CONSTRUCTORS(CooMatrix);

    ~CooMatrix() noexcept override = default;

    [[using gnu: flatten, leaf]] CooMatrix(const LilMatrix<T, Index>& lil_matrix) noexcept
        : nrows_{lil_matrix.nrows_}, ncols_{lil_matrix.ncols_} {
        const std::size_t values_size = lil_matrix.nnzs();
        offsets_map_.reserve(values_size);
        values_.reserve(values_size);
        row_indices_.reserve(values_size);
        col_indices_.reserve(values_size);
        for (const auto& [row, offsets_map] : lil_matrix.row_offsets_map_) {
            for (const auto& [col, offset] : offsets_map) {
                offsets_map_.emplace(IndexPair{row, col}, values_.size());
                values_.push_back(lil_matrix.row_values_map_.at(row)[offset]);
                row_indices_.push_back(row);
                col_indices_.push_back(col);
            }
        }
    }

    [[using gnu: pure, always_inline]]
    std::size_t nrows() const noexcept {
        return nrows_;
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
        if (nrows < nrows_ || ncols < ncols_) {
            const std::size_t values_size = values_.size();
            std::size_t j{0};
            for (std::size_t i{0}; i < values_size; ++i) {
                if (row_indices_[i] >= static_cast<Index>(nrows) ||
                    col_indices_[i] >= static_cast<Index>(ncols)) {
                    offsets_map_.erase(IndexPair{row_indices_[i], col_indices_[i]});
                    continue;
                }
                if (j != i) {
                    offsets_map_[IndexPair{row_indices_[i], col_indices_[i]}] = j;
                    values_[j] = values_[i];
                    row_indices_[j] = row_indices_[i];
                    col_indices_[j] = col_indices_[i];
                }
                ++j;
            }
            values_.resize(j);
            row_indices_.resize(j);
            col_indices_.resize(j);
            values_.shrink_to_fit();
            row_indices_.shrink_to_fit();
            col_indices_.shrink_to_fit();
        }
        nrows_ = nrows;
        ncols_ = ncols;
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    T coeff(Index row, Index col) const noexcept {
        if (row >= static_cast<Index>(nrows_) || col >= static_cast<Index>(ncols_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, nrows_, ncols_
            );
            return 0.0;
        }
        IndexPair index_pair{row, col};
        if (const auto search = offsets_map_.find(index_pair); search != offsets_map_.end()) {
            return values_[search->second];
        }
        return 0.0;
    }

    [[using gnu: always_inline, hot]]
    void updateCoeff(Index row, Index col, T coeff) noexcept {
        if (row >= static_cast<Index>(nrows_) || col >= static_cast<Index>(ncols_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, nrows_, ncols_
            );
            return;
        }
        IndexPair index_pair{row, col};
        if (auto search = offsets_map_.find(index_pair); search != offsets_map_.end()) {
            const std::size_t offset = search->second;
            values_[offset] = coeff;
        } else {
            if (coeff != 0.0) {
                offsets_map_[index_pair] = values_.size();
                values_.push_back(coeff);
                row_indices_.push_back(row);
                col_indices_.push_back(col);
            }
        }
        return;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& values() const noexcept {
        return values_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<Index>& rowIndices() const noexcept {
        return row_indices_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<Index>& colIndices() const noexcept {
        return col_indices_;
    }

    [[using gnu: always_inline]]
    void clear() noexcept {
        offsets_map_.clear();
        values_.clear();
        row_indices_.clear();
        col_indices_.clear();
        return;
    }

  private:
    struct [[nodiscard]] IndexPair final {
        Index row;
        Index col;

        [[using gnu: pure, always_inline]]
        bool
        operator==(IndexPair index_pair) const noexcept {
            return row == index_pair.row && col == index_pair.col;
        }

        template <typename Archive>
        [[using gnu: always_inline]]
        constexpr void serialize(Archive& ar, const unsigned int version) noexcept {
            ar& row;
            ar& col;
            return;
        }
    };

    struct [[nodiscard]] IndexPairHash final {
        [[using gnu: pure, always_inline]]
        std::size_t
        operator()(typename CooMatrix<T, Index>::IndexPair index_pair) const noexcept {
            return std::hash<Index>()(index_pair.row) ^ std::hash<Index>()(index_pair.col);
        }
    };

    struct [[nodiscard]] IndexPairRowMajorCompare final {
        [[using gnu: pure, always_inline, leaf]]
        bool
        operator()(const IndexPair& lhs, const IndexPair& rhs) const noexcept {
            if (lhs.row < rhs.row) {
                return true;
            } else if (lhs.row > rhs.row) {
                return false;
            }
            if (lhs.col < rhs.col) {
                return true;
            } else if (lhs.col > rhs.col) {
                return false;
            }
            return false;
        }
    };

    struct [[nodiscard]] IndexPairColumnMajorCompare final {
        [[using gnu: pure, always_inline, leaf]]
        bool
        operator()(const IndexPair& lhs, const IndexPair& rhs) const noexcept {
            if (lhs.col < rhs.col) {
                return true;
            } else if (lhs.col > rhs.col) {
                return false;
            }
            if (lhs.row < rhs.row) {
                return true;
            } else if (lhs.row > rhs.row) {
                return false;
            }
            return false;
        }
    };

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& nrows_;
        ar& ncols_;
        ar& offsets_map_;
        ar& values_;
        ar& row_indices_;
        ar& col_indices_;
        return;
    }

    std::size_t nrows_;
    std::size_t ncols_;
    std::unordered_map<
        typename CooMatrix<T, Index>::IndexPair, std::size_t,
        typename CooMatrix<T, Index>::IndexPairHash>
        offsets_map_;
    std::vector<T> values_;
    std::vector<Index> row_indices_;
    std::vector<Index> col_indices_;
};

} // namespace math
} // namespace tiny_pnc
