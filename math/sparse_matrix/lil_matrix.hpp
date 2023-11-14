/**
 * @file lil_matrix.hpp
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
class CooMatrix;

template <typename T, typename Index>
class CscMatrix;

template <typename T = double, typename Index = int>
class [[nodiscard]] LilMatrix final : public SparseMatrix<LilMatrix<T, Index>, T, Index> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_integral_v<Index>, "The loaded type must be a floating-point type.");

    friend class boost::serialization::access;
    friend class CooMatrix<T, Index>;

  public:
    [[using gnu: always_inline]] LilMatrix() noexcept
        : nrows_{}, ncols_{}, nnzs_{0}, row_offsets_map_{}, row_values_map_{}, row_indices_map_{} {}

    [[using gnu: always_inline]] explicit LilMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : nrows_{nrows}, ncols_{ncols}, nnzs_{0}, row_offsets_map_{}, row_values_map_{},
          row_indices_map_{} {}

    ENABLE_COPY_AND_MOVE(LilMatrix);

    ~LilMatrix() noexcept override = default;

    [[using gnu: flatten, leaf]] LilMatrix(const CooMatrix<T, Index>& coo_matrix) noexcept
        : nrows_{coo_matrix.nrows_}, ncols_{coo_matrix.ncols_}, nnzs_{coo_matrix.values_.size()} {
        for (const auto& [index_pair, offset] : coo_matrix.offsets_map_) {
            row_offsets_map_[index_pair.row].emplace(
                index_pair.col, row_values_map_[index_pair.row].size()
            );
            row_values_map_[index_pair.row].emplace_back(coo_matrix.values_[offset]);
            row_indices_map_[index_pair.row].emplace_back(index_pair.col);
        }
        for (auto& [row, values] : row_values_map_) {
            values.shrink_to_fit();
        }
        for (auto& [row, indices] : row_indices_map_) {
            indices.shrink_to_fit();
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
        return nnzs_;
    }

    [[using gnu: flatten, leaf]]
    void resize(std::size_t nrows, std::size_t ncols) noexcept {
        if (nrows < nrows_ || ncols < ncols_) {
            std::vector<Index> remove_rows;
            for (const auto& [row, row_offsets_map] : row_offsets_map_) {
                if (row >= static_cast<Index>(nrows)) {
                    remove_rows.push_back(row);
                }
            }
            for (const Index row : remove_rows) {
                row_offsets_map_.erase(row);
                nnzs_ -= row_values_map_[row].size();
                row_values_map_.erase(row);
                row_indices_map_.erase(row);
            }
            remove_rows.clear();
            for (auto& [row, indices] : row_indices_map_) {
                nnzs_ -= indices.size();
                std::size_t j{0};
                for (std::size_t i{0}; i < indices.size(); ++i) {
                    if (indices[i] >= static_cast<Index>(ncols)) {
                        row_offsets_map_[row].erase(indices[i]);
                        continue;
                    }
                    if (j != i) {
                        row_offsets_map_[row][indices[i]] = j;
                        row_values_map_[row][j] = row_values_map_[row][i];
                        indices[j] = indices[i];
                    }
                    ++j;
                }
                if (j != 0) {
                    row_values_map_[row].resize(j);
                    row_indices_map_[row].resize(j);
                    indices.shrink_to_fit();
                    indices.shrink_to_fit();
                    nnzs_ += j;
                } else {
                    remove_rows.push_back(row);
                }
            }
            for (const Index row : remove_rows) {
                row_offsets_map_.erase(row);
                row_values_map_.erase(row);
                row_indices_map_.erase(row);
            }
            remove_rows.clear();
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
        if (const auto search_row = row_offsets_map_.find(row);
            search_row != row_offsets_map_.end()) {
            if (const auto search_col = search_row->second.find(col);
                search_col != search_row->second.end()) {
                const std::size_t offset = search_col->second;
                return row_values_map_.at(row)[offset];
            }
        }
        return 0.0;
    }

    [[using gnu: hot]]
    void updateRow(Index row, std::vector<T> values, std::vector<Index> indices) noexcept {
        if (row >= static_cast<Index>(nrows_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" should be less than \"nrows\": row = {0:d} "
                "while nrows = {1:d}.",
                row, nrows_
            );
            return;
        }
        if (values.size() != indices.size()) {
            const std::size_t values_size{std::min(values.size(), indices.size())};
            values.resize(values_size);
            indices.resize(values_size);
        }
        if (auto search = row_offsets_map_.find(row); search != row_offsets_map_.end()) {
            row_offsets_map_.erase(search);
        }
        if (auto search = row_values_map_.find(row); search != row_values_map_.end()) {
            nnzs_ -= search->second.size();
            row_values_map_.erase(search);
        }
        if (auto search = row_indices_map_.find(row); search != row_indices_map_.end()) {
            row_indices_map_.erase(search);
        }
        if (values.size() == 0) {
            return;
        }
        const std::size_t values_size = values.size();
        std::unordered_map<Index, std::size_t> indices_map;
        indices_map.reserve(values_size);
        std::size_t j = 0;
        for (std::size_t i = 0; i < values_size; ++i) {
            if (indices[i] >= static_cast<Index>(ncols_)) {
                spdlog::warn(
                    "Out of range issue detected! \"col\" should be less than \"ncols\": col = "
                    "{0:d} while ncols = {1:d}.",
                    indices[i], ncols_
                );
                continue;
            }
            if (values[i] == 0.0) {
                continue;
            }
            indices_map.emplace(indices[i], j);
            if (j != i) {
                values[j] = values[i];
                indices[j] = indices[i];
            }
            ++j;
        }
        if (j != 0) {
            values.resize(j);
            indices.resize(j);
            values.shrink_to_fit();
            indices.shrink_to_fit();
            row_offsets_map_.emplace(row, std::move(indices_map));
            nnzs_ += values.size();
            row_values_map_.emplace(row, std::move(values));
            row_indices_map_.emplace(row, std::move(indices));
        }
        return;
    }

    [[using gnu: hot]]
    void updateRow(Index row, const std::unordered_map<Index, T>& values_map) noexcept {
        if (row >= static_cast<Index>(nrows_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" should be less than \"nrows\": row = {0:d} "
                "while nrows = {1:d}.",
                row, nrows_
            );
            return;
        }
        if (auto search = row_offsets_map_.find(row); search != row_offsets_map_.end()) {
            row_offsets_map_.erase(search);
        }
        if (auto search = row_values_map_.find(row); search != row_values_map_.end()) {
            nnzs_ -= search->second.size();
            row_values_map_.erase(search);
        }
        if (auto search = row_indices_map_.find(row); search != row_indices_map_.end()) {
            row_indices_map_.erase(search);
        }
        if (values_map.size() == 0) {
            return;
        }
        std::unordered_map<Index, std::size_t> indices_map;
        indices_map.reserve(values_map.size());
        std::vector<T> values;
        values.reserve(values_map.size());
        std::vector<Index> indices;
        indices.reserve(values_map.size());
        for (const auto& [col, value] : values_map) {
            if (col >= static_cast<Index>(ncols_)) {
                spdlog::warn(
                    "Out of range issue detected! \"col\" should be less than \"ncols\": col = "
                    "{0:d} while ncols = {1:d}.",
                    col, ncols_
                );
                continue;
            }
            if (value == 0.0) {
                continue;
            }
            indices_map.emplace(col, values.size());
            values.push_back(value);
            indices.push_back(col);
        }
        if (!values.empty()) {
            values.shrink_to_fit();
            indices.shrink_to_fit();
            row_offsets_map_.emplace(row, std::move(indices_map));
            nnzs_ += values.size();
            row_values_map_.emplace(row, std::move(values));
            row_indices_map_.emplace(row, std::move(indices));
        }
        return;
    }

    [[using gnu: pure, always_inline]]
    std::vector<T> rowValues(Index row) const noexcept {
        if (row >= static_cast<Index>(nrows_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" should be less than \"nrows\": row = {0:d} "
                "while nrows = {1:d}.",
                row, nrows_
            );
            return {};
        }
        if (const auto search = row_offsets_map_.find(row); search != row_offsets_map_.end()) {
            return row_values_map_.at(row);
        }
        return {};
    }

    [[using gnu: pure, always_inline]]
    std::vector<Index> rowIndices(Index row) const noexcept {
        if (row >= static_cast<Index>(nrows_)) {
            spdlog::warn(
                "Out of range issue detected! \"row\" should be less than \"nrows\": row = {0:d} "
                "while nrows = {1:d}.",
                row, nrows_
            );
            return {};
        }
        if (const auto search = row_offsets_map_.find(row); search != row_offsets_map_.end()) {
            return row_indices_map_.at(row);
        }
        return {};
    }

    [[using gnu: always_inline]]
    void clear() noexcept {
        nnzs_ = 0;
        row_offsets_map_.clear();
        row_values_map_.clear();
        row_indices_map_.clear();
        return;
    }

  private:
    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& nrows_;
        ar& ncols_;
        ar& nnzs_;
        ar& row_offsets_map_;
        ar& row_values_map_;
        ar& row_indices_map_;
        return;
    }

    std::size_t nrows_;
    std::size_t ncols_;
    std::size_t nnzs_;
    std::unordered_map<Index, std::unordered_map<Index, std::size_t>> row_offsets_map_;
    std::unordered_map<Index, std::vector<T>> row_values_map_;
    std::unordered_map<Index, std::vector<Index>> row_indices_map_;
};

} // namespace math
} // namespace tiny_pnc
