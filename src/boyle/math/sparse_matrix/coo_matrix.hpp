/**
 * @file coo_matrix.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-20
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <unordered_map>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common//utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse_matrix/dok_matrix.hpp"
#include "boyle/math/sparse_matrix/index_pair.hpp"

namespace boyle::math {

template <GeneralArithmetic Scalar = double, std::integral Index = int>
class CooMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;

    ENABLE_IMPLICIT_CONSTRUCTORS(CooMatrix);
    ~CooMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit CooMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

    [[using gnu: always_inline]]
    CooMatrix(const DokMatrix<Scalar, Index>& dok_matrix) noexcept
        : m_nrows{dok_matrix.nrows()}, m_ncols{dok_matrix.ncols()} {
        reserve(dok_matrix.nnzs());
        for (const auto& [index_pair, value] : dok_matrix.dictionary()) {
            if (value == Scalar{0.0}) {
                continue;
            }
            const auto& [row, col] = index_pair;
            m_offsets_map.emplace(index_pair, m_values.size());
            m_values.push_back(value);
            m_row_indices.push_back(row);
            m_col_indices.push_back(col);
        }
    }

    [[using gnu: flatten, leaf]] operator DokMatrix<Scalar, Index>() const noexcept {
        DokMatrix dok_matrix{m_nrows, m_ncols};
        const std::size_t nnzs{m_values.size()};
        dok_matrix.reserve(nnzs);
        for (std::size_t i{0}; i < nnzs; ++i) {
            dok_matrix.updateCoeff(m_row_indices[i], m_col_indices[i], m_values[i]);
        }
        return dok_matrix;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nrows() const noexcept -> std::size_t {
        return m_nrows;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto ncols() const noexcept -> std::size_t {
        return m_ncols;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nnzs() const noexcept -> std::size_t {
        return m_values.size();
    }

    [[using gnu: flatten, leaf]]
    auto resize(std::size_t nrows, std::size_t ncols) noexcept -> void {
        if (nrows < m_nrows || ncols < m_ncols) {
            const std::size_t nnzs{m_values.size()};
            std::size_t j{0};
            for (std::size_t i{0}; i < nnzs; ++i) {
                if (m_row_indices[i] >= static_cast<Index>(nrows) ||
                    m_col_indices[i] >= static_cast<Index>(ncols)) {
                    m_offsets_map.erase({m_row_indices[i], m_col_indices[i]});
                    continue;
                }
                m_offsets_map[{m_row_indices[i], m_col_indices[i]}] = j;
                m_values[j] = m_values[i];
                m_row_indices[j] = m_row_indices[i];
                m_col_indices[j] = m_col_indices[i];
                ++j;
            }
            m_values.resize(j);
            m_row_indices.resize(j);
            m_col_indices.resize(j);
            m_values.shrink_to_fit();
            m_row_indices.shrink_to_fit();
            m_col_indices.shrink_to_fit();
        }
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto reserve(std::size_t capacity) noexcept -> void {
        m_offsets_map.reserve(capacity);
        m_values.reserve(capacity);
        m_row_indices.reserve(capacity);
        m_col_indices.reserve(capacity);
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_offsets_map.clear();
        m_values.clear();
        m_row_indices.clear();
        m_col_indices.clear();
        return;
    }

    [[using gnu: flatten, leaf]]
    auto compress() noexcept -> void {
        const auto count = std::erase_if(m_offsets_map, [this](const auto& item) -> bool {
            return m_values[item.second] == Scalar{0.0};
        });
        if (count == 0) {
            return;
        }
        std::size_t j{0};
        for (const auto& [index_pair, offset] : m_offsets_map) {
            const auto& [row, col] = index_pair;
            m_values[j] = m_values[offset];
            m_row_indices[j] = row;
            m_col_indices[j] = col;
            ++j;
        }
        m_values.resize(j);
        m_row_indices.resize(j);
        m_col_indices.resize(j);
        m_values.shrink_to_fit();
        m_row_indices.shrink_to_fit();
        m_col_indices.shrink_to_fit();
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    auto coeff(Index row, Index col) const noexcept -> Scalar {
        if (row >= static_cast<Index>(m_nrows) || col >= static_cast<Index>(m_ncols)) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, m_nrows, m_ncols
            );
            return 0.0;
        }
        if (const auto search = m_offsets_map.find({row, col}); search != m_offsets_map.cend()) {
            return m_values[search->second];
        }
        return 0.0;
    }

    [[using gnu: always_inline, hot]]
    auto updateCoeff(Index row, Index col, Scalar value) noexcept -> void {
        if (row >= static_cast<Index>(m_nrows) || col >= static_cast<Index>(m_ncols)) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, m_nrows, m_ncols
            );
            return;
        }
        if (const auto search = m_offsets_map.find({row, col}); search != m_offsets_map.cend()) {
            m_values[search->second] = value;
        } else {
            if (value != 0.0) {
                m_offsets_map.emplace(IndexPair<Index>{row, col}, m_values.size());
                m_values.push_back(value);
                m_row_indices.push_back(row);
                m_col_indices.push_back(col);
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(Index row, Index col) const noexcept -> Scalar {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto values() const noexcept -> const std::vector<Scalar>& {
        return m_values;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto rowIndices() const noexcept -> const std::vector<Index>& {
        return m_row_indices;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto colIndices() const noexcept -> const std::vector<Index>& {
        return m_col_indices;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_nrows;
        archive & m_ncols;
        archive & m_offsets_map;
        archive & m_values;
        archive & m_row_indices;
        archive & m_col_indices;
        return;
    }

    std::size_t m_nrows{0};
    std::size_t m_ncols{0};
    std::unordered_map<IndexPair<Index>, std::size_t, IndexPairHash<Index>> m_offsets_map{};
    std::vector<Scalar> m_values{};
    std::vector<Index> m_row_indices{};
    std::vector<Index> m_col_indices{};
};

} // namespace boyle::math
