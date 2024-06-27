/**
 * @file csr_matrix.hpp
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
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse_matrix/dok_matrix.hpp"
#include "boyle/math/sparse_matrix/index_pair.hpp"

namespace boyle::math {

template <GeneralArithmetic Scalar = double, std::integral Index = int>
class [[nodiscard]] CsrMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;

    ENABLE_IMPLICIT_CONSTRUCTORS(CsrMatrix);
    ~CsrMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit CsrMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : m_ncols{ncols}, m_outer_indices(nrows + 1, 0) {}

    [[using gnu: flatten, leaf]]
    CsrMatrix(const DokMatrix<Scalar, Index>& dok_matrix) noexcept
        : m_ncols{dok_matrix.ncols()} {
        const std::size_t nnzs = dok_matrix.nnzs();
        reserve(nnzs);
        m_outer_indices.reserve(dok_matrix.nrows() + 1);
        m_outer_indices.push_back(0);
        std::vector<std::pair<IndexPair<Index>, Scalar>> row_major_dictionary{
            dok_matrix.dictionary().cbegin(), dok_matrix.dictionary().cend()
        };
        std::sort(
            row_major_dictionary.begin(), row_major_dictionary.end(),
            [](const std::pair<IndexPair<Index>, Scalar>& a,
               const std::pair<IndexPair<Index>, Scalar>& b) noexcept {
                constexpr IndexPairRowMajorCompare<Index> Compare{};
                return Compare(a.first, b.first);
            }
        );
        for (const auto& [index_pair, value] : row_major_dictionary) {
            const auto& [row, col] = index_pair;
            m_values.push_back(value);
            m_inner_indices.push_back(col);
            const std::size_t diff_row = row - m_outer_indices.size() + 1;
            for (std::size_t i{0}; i < diff_row; ++i) {
                m_outer_indices.push_back(m_values.size() - 1);
            }
        }
        for (std::size_t i{m_outer_indices.size()}; i < dok_matrix.nrows() + 1; ++i) {
            m_outer_indices.push_back(nnzs);
        }
    }

    [[using gnu: flatten, leaf]] operator DokMatrix<Scalar, Index>() const noexcept {
        const std::size_t nrows{m_outer_indices.size() - 1};
        const std::size_t ncols{m_ncols};
        DokMatrix<Scalar, Index> dok_matrix{nrows, ncols};
        for (Index row{0}; row < static_cast<Index>(nrows); ++row) {
            for (Index offset{m_outer_indices[row]}; offset < m_outer_indices[row + 1]; ++offset) {
                const Index col{m_inner_indices[offset]};
                dok_matrix.updateCoeff(row, col, m_values[offset]);
            }
        }
        return dok_matrix;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nrows() const noexcept -> std::size_t {
        return m_outer_indices.size() - 1;
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
        if (nrows >= m_outer_indices.size() - 1) {
            for (std::size_t i{m_outer_indices.size() - 1}; i < nrows; ++i) {
                m_outer_indices.push_back(m_outer_indices.back());
            }
        } else {
            m_outer_indices.resize(nrows + 1);
            m_values.resize(m_outer_indices.back());
            m_inner_indices.resize(m_outer_indices.back());
        }
        const std::size_t nnzs{m_values.size()};
        if (ncols < m_ncols) {
            std::size_t j{0}, k{0};
            m_ncols = ncols;
            for (std::size_t i{0}; i < nnzs; ++i) {
                if (m_inner_indices[i] >= static_cast<Index>(m_ncols)) {
                    continue;
                }
                if (j != i) {
                    m_values[j] = m_values[i];
                    m_inner_indices[j] = m_inner_indices[i];
                    if (i == static_cast<std::size_t>(m_outer_indices[k])) {
                        m_outer_indices[k] = j;
                    }
                }
                ++j;
                if (i == static_cast<std::size_t>(m_outer_indices[k])) {
                    ++k;
                }
            }
            m_outer_indices.back() = j;
            m_values.resize(j);
            m_inner_indices.resize(j);
        }
        return;
    }

    [[using gnu: always_inline]]
    auto reserve(std::size_t capacity) noexcept -> void {
        m_values.reserve(capacity);
        m_inner_indices.reserve(capacity);
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_values.clear();
        m_inner_indices.clear();
        std::fill(m_outer_indices.begin(), m_outer_indices.end(), 0);
        return;
    }

    [[using gnu: flatten, leaf]]
    auto compress() noexcept -> void {
        const std::size_t nnzs{m_values.size()};
        std::size_t j{0}, k{0};
        for (std::size_t i{0}; i < nnzs; ++i) {
            if (m_values[i] == Scalar{0.0}) {
                continue;
            }
            if (j != i) {
                m_values[j] = m_values[i];
                m_inner_indices[j] = m_inner_indices[i];
                if (i == static_cast<std::size_t>(m_outer_indices[k])) {
                    m_outer_indices[k] = j;
                }
            }
            ++j;
            if (i == static_cast<std::size_t>(m_outer_indices[k])) {
                ++k;
            }
        }
        m_outer_indices.back() = j;
        m_values.resize(j);
        m_inner_indices.resize(j);
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    auto coeff(Index row, Index col) const noexcept -> Scalar {
        if (row >= static_cast<Index>(m_outer_indices.size() - 1) ||
            col >= static_cast<Index>(m_ncols)) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, m_outer_indices.size() - 1, m_ncols
            );
            return 0.0;
        }
        if (const auto search = std::find(
                m_inner_indices.cbegin() + m_outer_indices[row],
                m_inner_indices.cbegin() + m_outer_indices[row + 1], col
            );
            search != m_inner_indices.cbegin() + m_outer_indices[row + 1]) {
            return m_values[search - m_inner_indices.cbegin()];
        }
        return 0.0;
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(Index row, Index col, Scalar value) noexcept -> void {
        if (row >= static_cast<Index>(m_outer_indices.size() - 1) ||
            col >= static_cast<Index>(m_ncols)) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! \"row\" and \"col\" should be less than \"nrows\" "
                "and \"ncols\" respectively: (row, col) = ({0:d}, {1:d}) while (nrows, ncols) = "
                "({2:d}, {3:d}).",
                row, col, m_outer_indices.size() - 1, m_ncols
            );
            return;
        }
        if (const auto search = std::find(
                m_inner_indices.cbegin() + m_outer_indices[row],
                m_inner_indices.cbegin() + m_outer_indices[row + 1], col
            );
            search != m_inner_indices.cbegin() + m_outer_indices[row + 1]) {
            m_values[search - m_inner_indices.cbegin()] = value;
        } else {
            if (value != 0.0) {
                const std::size_t inner_offset =
                    std::ranges::upper_bound(
                        m_inner_indices.cbegin() + m_outer_indices[row],
                        m_inner_indices.cbegin() + m_outer_indices[row + 1], col
                    ) -
                    m_inner_indices.cbegin();
                m_values.insert(m_values.begin() + inner_offset, value);
                m_inner_indices.insert(m_inner_indices.begin() + inner_offset, col);
                for (auto outer_it = m_outer_indices.begin() + row + 1;
                     outer_it != m_outer_indices.end(); ++outer_it) {
                    ++(*outer_it);
                }
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
    auto innerIndices() const noexcept -> const std::vector<Index>& {
        return m_inner_indices;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto outerIndices() const noexcept -> const std::vector<Index>& {
        return m_outer_indices;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_ncols;
        archive & m_values;
        archive & m_inner_indices;
        archive & m_outer_indices;
        return;
    }

    std::size_t m_ncols{0};
    std::vector<Scalar> m_values{};
    std::vector<Index> m_inner_indices{};
    std::vector<Index> m_outer_indices{};
};

} // namespace boyle::math
