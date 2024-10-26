/**
 * @file lil_matrix.hpp
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

#include <algorithm>
#include <concepts>
#include <utility>

#include "boost/serialization/access.hpp"
#include "boost/unordered/unordered_flat_map.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/dok_matrix.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar = double, std::integral Index = int>
class [[nodiscard]] LilMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    LilMatrix() noexcept = default;
    LilMatrix(const LilMatrix& other) noexcept = default;
    auto operator=(const LilMatrix& other) noexcept -> LilMatrix& = default;
    LilMatrix(LilMatrix&& other) noexcept = default;
    auto operator=(LilMatrix&& other) noexcept -> LilMatrix& = default;
    ~LilMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit LilMatrix(size_type nrows, size_type ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

    [[using gnu: always_inline]]
    LilMatrix(const DokMatrix<value_type, index_type>& dok_matrix) noexcept
        : m_nrows{dok_matrix.nrows()}, m_ncols{dok_matrix.ncols()}, m_nnzs{dok_matrix.nnzs()} {
        for (const auto& [index_pair, value] : dok_matrix.dictionary()) {
            m_row_dictionaries[index_pair.row].emplace(index_pair.col, value);
        }
    }

    [[using gnu: ]] operator DokMatrix<value_type, index_type>() const noexcept {
        DokMatrix<value_type, index_type> dok_matrix{m_nrows, m_ncols};
        dok_matrix.reserve(m_nnzs);
        for (const auto& [row, row_dictionary] : m_row_dictionaries) {
            for (const auto& [col, value] : row_dictionary) {
                dok_matrix.updateCoeff(row, col, value);
            }
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

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto nnzs() const noexcept -> size_type {
        return m_nnzs;
    }

    [[using gnu: ]]
    auto resize(size_type nrows, size_type ncols) noexcept -> void {
        boost::unordered::erase_if(
            m_row_dictionaries,
            [nrows](const auto& row_dictionary) noexcept -> bool {
                return row_dictionary.first >= static_cast<index_type>(nrows);
            }
        );
        for (auto& [row, row_dictionary] : m_row_dictionaries) {
            boost::unordered::erase_if(row_dictionary, [ncols](const auto& item) noexcept -> bool {
                return item.first >= static_cast<index_type>(ncols);
            });
        }
        boost::unordered::erase_if(
            m_row_dictionaries,
            [](const auto& row_dictionary) noexcept -> bool {
                return row_dictionary.second.empty();
            }
        );
        m_nnzs = 0;
        for (const auto& [row, row_dictionary] : m_row_dictionaries) {
            m_nnzs += row_dictionary.size();
        }
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto reserve(size_type capacity) noexcept -> void {
        m_row_dictionaries.reserve(capacity);
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_nnzs = 0;
        m_row_dictionaries.clear();
        return;
    }

    [[using gnu: always_inline]]
    auto compress() noexcept -> void {
        return;
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto coeff(size_type row, size_type col) const noexcept -> value_type {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return static_cast<value_type>(0.0);
        }
        if (const auto search_row = m_row_dictionaries.find(static_cast<index_type>(row));
            search_row != m_row_dictionaries.cend()) {
            if (const auto search_col = search_row->second.find(static_cast<index_type>(col));
                search_col != search_row->second.cend()) {
                return search_col->second;
            }
        }
        return static_cast<value_type>(0.0);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(size_type row, size_type col, const_reference value) noexcept -> void {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return;
        }
        if (auto search_row = m_row_dictionaries.find(static_cast<index_type>(row));
            search_row != m_row_dictionaries.end()) {
            if (auto search_col = search_row->second.find(static_cast<index_type>(col));
                search_col != search_row->second.end()) {
                if (value != value_type(0.0)) {
                    search_col->second = value;
                } else {
                    m_nnzs -= 1;
                    search_row->second.erase(search_col);
                    if (search_row->second.empty()) {
                        m_row_dictionaries.erase(search_row);
                    }
                }
            } else {
                if (value != value_type(0.0)) {
                    m_nnzs += 1;
                    search_row->second.emplace(static_cast<index_type>(col), value);
                }
            }
        } else {
            if (value != value_type(0.0)) {
                m_nnzs += 1;
                m_row_dictionaries.emplace(
                    static_cast<index_type>(row),
                    boost::unordered_flat_map<int, double>{{static_cast<index_type>(col), value}}
                );
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type row, size_type col) const noexcept -> value_type {
        return coeff(row, col);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateRow(
        size_type row, boost::unordered_flat_map<index_type, value_type> row_dictionary
    ) noexcept -> void {
        if (row >= m_nrows) [[unlikely]] {
            return;
        }
        if (auto search = m_row_dictionaries.find(static_cast<index_type>(row));
            search != m_row_dictionaries.end()) {
            m_nnzs -= search->second.size();
            m_row_dictionaries.erase(search);
        }
        boost::unordered::erase_if(row_dictionary, [this](const auto& item) noexcept -> bool {
            return item.first >= static_cast<index_type>(m_ncols) || item.second == value_type(0.0);
        });
        if (row_dictionary.empty()) {
            return;
        }
        m_nnzs += row_dictionary.size();
        m_row_dictionaries.emplace(static_cast<index_type>(row), std::move(row_dictionary));
        return;
    }

    [[using gnu: always_inline, leaf]]
    auto row_dictionaries() const noexcept -> const
        boost::unordered_flat_map<index_type, boost::unordered_flat_map<index_type, value_type>>& {
        return m_row_dictionaries;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_nrows;
        archive & m_ncols;
        archive & m_nnzs;
        archive & m_row_dictionaries;
        return;
    }

    size_type m_nrows{0};
    size_type m_ncols{0};
    size_type m_nnzs{0};
    boost::unordered_flat_map<index_type, boost::unordered_flat_map<index_type, value_type>>
        m_row_dictionaries{};
};

} // namespace boyle::math
