/**
 * @file csc_matrix.hpp
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
#include <span>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/sparse_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, std::integral Index, typename Alloc>
class [[nodiscard]] CscMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = typename SparseTraits<CscMatrix>::value_type;
    using index_type = typename SparseTraits<CscMatrix>::index_type;
    using reference = typename SparseTraits<CscMatrix>::reference;
    using const_reference = typename SparseTraits<CscMatrix>::const_reference;
    using size_type = typename SparseTraits<CscMatrix>::size_type;
    using difference_type = typename SparseTraits<CscMatrix>::difference_type;
    using allocator_type = typename SparseTraits<CscMatrix>::allocator_type;

    [[using gnu: always_inline]]
    explicit CscMatrix(const allocator_type& alloc = {}) noexcept
        : m_values(alloc),
          m_inner_indices(
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>(
                  alloc
              )
          ),
          m_outer_indices(
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>(
                  alloc
              )
          ) {}
    CscMatrix(const CscMatrix& other) noexcept = default;
    auto operator=(const CscMatrix& other) noexcept -> CscMatrix& = default;
    CscMatrix(CscMatrix&& other) noexcept = default;
    auto operator=(CscMatrix&& other) noexcept -> CscMatrix& = default;
    ~CscMatrix() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_values.get_allocator();
    }

    [[using gnu: always_inline]]
    explicit CscMatrix(size_type nrows, size_type ncols, const allocator_type& alloc = {}) noexcept
        : m_nrows{nrows}, m_values(alloc),
          m_inner_indices(
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>(
                  alloc
              )
          ),
          m_outer_indices(
              ncols + 1, 0,
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>(
                  alloc
              )
          ) {}

    [[using gnu: ]]
    CscMatrix(const DokMatrix<value_type, index_type, allocator_type>& dok_matrix) noexcept
        : m_nrows{dok_matrix.nrows()}, m_values(dok_matrix.get_allocator()),
          m_inner_indices(
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>{
                  dok_matrix.get_allocator()
              }
          ),
          m_outer_indices(
              typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>{
                  dok_matrix.get_allocator()
              }
          ) {
        const size_type nnzs{dok_matrix.nnzs()};
        reserve(nnzs);
        m_outer_indices.reserve(dok_matrix.ncols() + 1);
        m_outer_indices.push_back(0);
        static std::array<std::pair<IndexPair<index_type>, value_type>, 2048>
            column_major_dictionary;
        auto [_, column_major_dictionary_end] = std::ranges::partial_sort_copy(
            dok_matrix.dictionary().cbegin(), dok_matrix.dictionary().cend(),
            column_major_dictionary.begin(), column_major_dictionary.end(),
            [](const std::pair<IndexPair<index_type>, value_type>& a,
               const std::pair<IndexPair<index_type>, value_type>& b) static constexpr noexcept
                -> bool {
                constexpr IndexPairColumnMajorCompare<index_type> Compare{};
                return Compare(a.first, b.first);
            }
        );
        // std::vector<
        //     std::pair<IndexPair<index_type>, value_type>,
        //     typename std::allocator_traits<allocator_type>::template rebind_alloc<
        //         std::pair<IndexPair<index_type>, value_type>>>
        //     column_major_dictionary{
        //         dok_matrix.dictionary().cbegin(), dok_matrix.dictionary().cend(),
        //         typename std::allocator_traits<allocator_type>::template rebind_alloc<
        //             std::pair<IndexPair<index_type>, value_type>>{dok_matrix.get_allocator()}
        //     };
        // std::ranges::sort(
        //     column_major_dictionary.begin(), column_major_dictionary.end(),
        //     [](const std::pair<IndexPair<index_type>, value_type>& a,
        //        const std::pair<IndexPair<index_type>, value_type>& b) static constexpr noexcept
        //        -> bool {
        //         constexpr IndexPairColumnMajorCompare<index_type> Compare{};
        //         return Compare(a.first, b.first);
        //     }
        // );
        for (const auto& [index_pair, value] :
             std::span{column_major_dictionary.begin(), column_major_dictionary_end}) {
            const auto& [row, col] = index_pair;
            m_values.push_back(value);
            m_inner_indices.push_back(row);
            const size_type diff_col = col - m_outer_indices.size() + 1;
            for (size_type i{0}; i < diff_col; ++i) {
                m_outer_indices.push_back(m_values.size() - 1);
            }
        }
        for (size_type i{m_outer_indices.size()}; i < dok_matrix.ncols() + 1; ++i) {
            m_outer_indices.push_back(nnzs);
        }
    }

    [[using gnu: ]] operator DokMatrix<value_type, index_type, allocator_type>() const noexcept {
        const size_type nrows{m_nrows};
        const size_type ncols{m_outer_indices.size() - 1};
        DokMatrix<value_type, index_type, allocator_type> dok_matrix{nrows, ncols, get_allocator()};
        for (size_type j{0}; j < ncols; ++j) {
            for (size_type offset(m_outer_indices[j]);
                 offset < static_cast<size_type>(m_outer_indices[j + 1]); ++offset) {
                const size_type i(m_inner_indices[offset]);
                dok_matrix.updateCoeff(i, j, m_values[offset]);
            }
        }
        return dok_matrix;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto nrows() const noexcept -> size_type {
        return m_nrows;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto ncols() const noexcept -> size_type {
        return m_outer_indices.size() - 1;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nnzs() const noexcept -> size_type {
        return m_values.size();
    }

    [[using gnu: ]]
    auto resize(size_type nrows, size_type ncols) noexcept -> void {
        if (ncols >= m_outer_indices.size() - 1) {
            for (size_type i{m_outer_indices.size() - 1}; i < ncols; ++i) {
                m_outer_indices.push_back(m_outer_indices.back());
            }
        } else {
            m_outer_indices.resize(ncols + 1);
            m_values.resize(m_outer_indices.back());
            m_inner_indices.resize(m_outer_indices.back());
        }
        const size_type nnzs{m_values.size()};
        if (nrows < m_nrows) {
            size_type j{0}, k{0};
            m_nrows = nrows;
            for (size_type i{0}; i < nnzs; ++i) {
                if (m_inner_indices[i] >= static_cast<index_type>(m_nrows)) {
                    continue;
                }
                if (j != i) {
                    m_values[j] = m_values[i];
                    m_inner_indices[j] = m_inner_indices[i];
                    if (i == static_cast<size_type>(m_outer_indices[k])) {
                        m_outer_indices[k] = j;
                    }
                }
                ++j;
                if (i == static_cast<size_type>(m_outer_indices[k])) {
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
    auto reserve(size_type capacity) noexcept -> void {
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

    [[using gnu: ]]
    auto compress() noexcept -> void {
        const size_type nnzs{m_values.size()};
        size_type j{0}, k{0};
        for (size_type i{0}; i < nnzs; ++i) {
            if (m_values[i] == value_type(0.0)) {
                continue;
            }
            if (j != i) {
                m_values[j] = m_values[i];
                m_inner_indices[j] = m_inner_indices[i];
                if (i == static_cast<size_type>(m_outer_indices[k])) {
                    m_outer_indices[k] = j;
                }
            }
            ++j;
            if (i == static_cast<size_type>(m_outer_indices[k])) {
                ++k;
            }
        }
        m_outer_indices.back() = j;
        m_values.resize(j);
        m_inner_indices.resize(j);
        return;
    }

    [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
    auto coeff(size_type row, size_type col) const noexcept -> value_type {
        if (row >= m_nrows || col + 1 >= m_outer_indices.size()) [[unlikely]] {
            return static_cast<value_type>(0.0);
        }
        if (const auto search = std::find(
                m_inner_indices.cbegin() + m_outer_indices[col],
                m_inner_indices.cbegin() + m_outer_indices[col + 1], static_cast<index_type>(row)
            );
            search != m_inner_indices.cbegin() + m_outer_indices[col + 1]) {
            return m_values[search - m_inner_indices.cbegin()];
        }
        return static_cast<value_type>(0.0);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(size_type row, size_type col, const_reference value) noexcept -> void {
        if (row >= m_nrows || col + 1 >= m_outer_indices.size()) [[unlikely]] {
            return;
        }
        if (const auto search = std::find(
                m_inner_indices.cbegin() + m_outer_indices[col],
                m_inner_indices.cbegin() + m_outer_indices[col + 1], static_cast<index_type>(row)
            );
            search != m_inner_indices.cbegin() + m_outer_indices[col + 1]) {
            m_values[search - m_inner_indices.cbegin()] = value;
        } else {
            const size_type inner_offset = std::ranges::upper_bound(
                                               m_inner_indices.cbegin() + m_outer_indices[col],
                                               m_inner_indices.cbegin() + m_outer_indices[col + 1],
                                               static_cast<index_type>(row)
                                           ) -
                                           m_inner_indices.cbegin();
            m_values.insert(m_values.begin() + inner_offset, value);
            m_inner_indices.insert(
                m_inner_indices.begin() + inner_offset, static_cast<index_type>(row)
            );
            for (auto outer_it = m_outer_indices.begin() + col + 1;
                 outer_it != m_outer_indices.end(); ++outer_it) {
                ++(*outer_it);
            }
        }
        return;
    }

    [[using gnu: flatten, leaf, hot]]
    auto eraseCoeff(size_type row, size_type col) noexcept -> void {
        if (row >= m_nrows || col + 1 >= m_outer_indices.size()) [[unlikely]] {
            return;
        }
        if (const auto search = std::find(
                m_inner_indices.cbegin() + m_outer_indices[col],
                m_inner_indices.cbegin() + m_outer_indices[col + 1], static_cast<index_type>(row)
            );
            search != m_inner_indices.cbegin() + m_outer_indices[col + 1]) {
            const size_type inner_offset = search - m_inner_indices.cbegin();
            m_values.erase(m_values.begin() + inner_offset);
            m_inner_indices.erase(m_inner_indices.begin() + inner_offset);
            for (auto outer_it = m_outer_indices.begin() + col + 1;
                 outer_it != m_outer_indices.end(); ++outer_it) {
                --(*outer_it);
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type row, size_type col) const noexcept -> value_type {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto values() const noexcept -> std::span<const value_type> {
        return m_values;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto innerIndices() const noexcept -> std::span<const index_type> {
        return m_inner_indices;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto outerIndices() const noexcept -> std::span<const index_type> {
        return m_outer_indices;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_nrows;
        archive & m_values;
        archive & m_inner_indices;
        archive & m_outer_indices;
        return;
    }

    size_type m_nrows{0};
    std::vector<value_type, allocator_type> m_values{};
    std::vector<
        index_type,
        typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>>
        m_inner_indices{};
    std::vector<
        index_type,
        typename std::allocator_traits<allocator_type>::template rebind_alloc<index_type>>
        m_outer_indices{};
};

} // namespace boyle::math
