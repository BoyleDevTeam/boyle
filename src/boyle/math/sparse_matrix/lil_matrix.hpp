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
#include "boyle/math/sparse_matrix/dok_matrix.hpp"

namespace boyle::math {

template <GeneralArithmetic Scalar = double, std::integral Index = int>
class [[nodiscard]] LilMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;

    LilMatrix() noexcept = default;
    LilMatrix(const LilMatrix& other) noexcept = default;
    auto operator=(const LilMatrix& other) noexcept -> LilMatrix& = default;
    LilMatrix(LilMatrix&& other) noexcept = default;
    auto operator=(LilMatrix&& other) noexcept -> LilMatrix& = default;
    ~LilMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit LilMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

    [[using gnu: always_inline]]
    LilMatrix(const DokMatrix<Scalar, Index>& dok_matrix) noexcept
        : m_nrows{dok_matrix.nrows()}, m_ncols{dok_matrix.ncols()}, m_nnzs{dok_matrix.nnzs()} {
        for (const auto& [index_pair, value] : dok_matrix.dictionary()) {
            m_row_dictionaries[index_pair.row].emplace(index_pair.col, value);
        }
    }

    [[using gnu: ]] operator DokMatrix<Scalar, Index>() const noexcept {
        DokMatrix dok_matrix{m_nrows, m_ncols};
        dok_matrix.reserve(m_nnzs);
        for (const auto& [row, row_dictionary] : m_row_dictionaries) {
            for (const auto& [col, value] : row_dictionary) {
                dok_matrix.updateCoeff(row, col, value);
            }
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
        return m_nnzs;
    }

    [[using gnu: ]]
    auto resize(std::size_t nrows, std::size_t ncols) noexcept -> void {
        boost::unordered::erase_if(
            m_row_dictionaries,
            [nrows](const auto& row_dictionary) noexcept -> bool {
                return row_dictionary.first >= static_cast<Index>(nrows);
            }
        );
        for (auto& [row, row_dictionary] : m_row_dictionaries) {
            boost::unordered::erase_if(row_dictionary, [ncols](const auto& item) noexcept -> bool {
                return item.first >= static_cast<Index>(ncols);
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
    auto reserve(std::size_t capacity) noexcept -> void {
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
    auto coeff(Index row, Index col) const noexcept -> Scalar {
        if (row >= static_cast<Index>(m_nrows) || col >= static_cast<Index>(m_ncols)) [[unlikely]] {
            return 0.0;
        }
        if (const auto search_row = m_row_dictionaries.find(row);
            search_row != m_row_dictionaries.cend()) {
            if (const auto search_col = search_row->second.find(col);
                search_col != search_row->second.cend()) {
                return search_col->second;
            }
        }
        return 0.0;
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(Index row, Index col, Scalar value) noexcept -> void {
        if (row >= static_cast<Index>(m_nrows) || col >= static_cast<Index>(m_ncols)) [[unlikely]] {
            return;
        }
        if (auto search_row = m_row_dictionaries.find(row);
            search_row != m_row_dictionaries.end()) {
            if (auto search_col = search_row->second.find(col);
                search_col != search_row->second.end()) {
                if (value != Scalar{0.0}) {
                    search_col->second = value;
                } else {
                    m_nnzs -= 1;
                    search_row->second.erase(search_col);
                    if (search_row->second.empty()) {
                        m_row_dictionaries.erase(search_row);
                    }
                }
            } else {
                if (value != Scalar{0.0}) {
                    m_nnzs += 1;
                    search_row->second.emplace(col, value);
                }
            }
        } else {
            if (value != Scalar{0.0}) {
                m_nnzs += 1;
                m_row_dictionaries.emplace(
                    row, boost::unordered_flat_map<int, double>{{col, value}}
                );
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(Index row, Index col) const noexcept -> Scalar {
        return coeff(row, col);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateRow(Index row, boost::unordered_flat_map<Index, Scalar> row_dictionary) noexcept
        -> void {
        if (row >= static_cast<Index>(m_nrows)) {
            return;
        }
        if (auto search = m_row_dictionaries.find(row); search != m_row_dictionaries.end()) {
            m_nnzs -= search->second.size();
            m_row_dictionaries.erase(search);
        }
        boost::unordered::erase_if(row_dictionary, [this](const auto& item) noexcept -> bool {
            return item.first >= static_cast<Index>(m_ncols) || item.second == Scalar{0.0};
        });
        if (row_dictionary.empty()) {
            return;
        }
        m_nnzs += row_dictionary.size();
        m_row_dictionaries.emplace(row, std::move(row_dictionary));
        return;
    }

    [[using gnu: always_inline, leaf]]
    auto row_dictionaries() const noexcept
        -> const boost::unordered_flat_map<Index, boost::unordered_flat_map<Index, Scalar>>& {
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

    std::size_t m_nrows{0};
    std::size_t m_ncols{0};
    std::size_t m_nnzs{0};
    boost::unordered_flat_map<Index, boost::unordered_flat_map<Index, Scalar>> m_row_dictionaries{};
};

} // namespace boyle::math
