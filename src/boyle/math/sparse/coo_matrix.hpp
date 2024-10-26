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
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/unordered/unordered_flat_map.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/dok_matrix.hpp"
#include "boyle/math/sparse/index_pair.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar = double, std::integral Index = int>
class CooMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    CooMatrix() noexcept = default;
    CooMatrix(const CooMatrix& other) noexcept = default;
    auto operator=(const CooMatrix& other) noexcept -> CooMatrix& = default;
    CooMatrix(CooMatrix&& other) noexcept = default;
    auto operator=(CooMatrix&& other) noexcept -> CooMatrix& = default;
    ~CooMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit CooMatrix(size_type nrows, size_type ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

    [[using gnu: ]]
    CooMatrix(const DokMatrix<value_type, index_type>& dok_matrix) noexcept
        : m_nrows{dok_matrix.nrows()}, m_ncols{dok_matrix.ncols()} {
        reserve(dok_matrix.nnzs());
        for (const auto& [index_pair, value] : dok_matrix.dictionary()) {
            if (value == static_cast<value_type>(0.0)) {
                continue;
            }
            const auto& [row, col] = index_pair;
            m_offsets_map.emplace(index_pair, m_values.size());
            m_values.push_back(value);
            m_row_indices.push_back(row);
            m_col_indices.push_back(col);
        }
    }

    [[using gnu: ]] operator DokMatrix<value_type, index_type>() const noexcept {
        DokMatrix<value_type, index_type> dok_matrix(m_nrows, m_ncols);
        const size_type nnzs{m_values.size()};
        dok_matrix.reserve(nnzs);
        for (size_type i{0}; i < nnzs; ++i) {
            dok_matrix.updateCoeff(m_row_indices[i], m_col_indices[i], m_values[i]);
        }
        return dok_matrix;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto nrows() const noexcept -> size_type {
        return m_nrows;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto ncols() const noexcept -> size_type {
        return m_ncols;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nnzs() const noexcept -> size_type {
        return m_values.size();
    }

    [[using gnu: ]]
    auto resize(size_type nrows, size_type ncols) noexcept -> void {
        if (nrows < m_nrows || ncols < m_ncols) {
            const size_type nnzs{m_values.size()};
            size_type j{0};
            for (size_type i{0}; i < nnzs; ++i) {
                if (m_row_indices[i] >= static_cast<index_type>(nrows) ||
                    m_col_indices[i] >= static_cast<index_type>(ncols)) {
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
    auto reserve(size_type capacity) noexcept -> void {
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

    [[using gnu: ]]
    auto compress() noexcept -> void {
        const auto count =
            boost::unordered::erase_if(m_offsets_map, [this](const auto& item) noexcept -> bool {
                return m_values[item.second] == static_cast<value_type>(0.0);
            });
        if (count == 0) {
            return;
        }
        size_type j{0};
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

    [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
    auto coeff(size_type row, size_type col) const noexcept -> value_type {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return static_cast<value_type>(0.0);
        }
        if (const auto search =
                m_offsets_map.find({static_cast<index_type>(row), static_cast<index_type>(col)});
            search != m_offsets_map.cend()) {
            return m_values[search->second];
        }
        return static_cast<value_type>(0.0);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(size_type row, size_type col, const_reference value) noexcept -> void {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return;
        }
        if (const auto search =
                m_offsets_map.find({static_cast<index_type>(row), static_cast<index_type>(col)});
            search != m_offsets_map.cend()) {
            m_values[search->second] = value;
        } else {
            if (value != static_cast<value_type>(0.0)) {
                m_offsets_map.emplace(
                    IndexPair<index_type>{
                        static_cast<index_type>(row), static_cast<index_type>(col)
                    },
                    m_values.size()
                );
                m_values.push_back(value);
                m_row_indices.push_back(row);
                m_col_indices.push_back(col);
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type row, size_type col) const noexcept -> value_type {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto values() const noexcept -> const std::vector<Scalar>& {
        return m_values;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto rowIndices() const noexcept -> const std::vector<Index>& {
        return m_row_indices;
    }

    [[using gnu: pure, always_inline, leaf]]
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

    size_type m_nrows{0};
    size_type m_ncols{0};
    boost::unordered_flat_map<
        IndexPair<index_type>, size_type, IndexPairHash<index_type>, IndexPairEqual<index_type>>
        m_offsets_map{};
    std::vector<value_type> m_values{};
    std::vector<index_type> m_row_indices{};
    std::vector<index_type> m_col_indices{};
};

} // namespace boyle::math
