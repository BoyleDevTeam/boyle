/**
 * @file dok_matrix.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-06-24
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>

#include "boost/serialization/access.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/sparse_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, std::integral Index, typename Alloc>
class DokMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = typename SparseTraits<DokMatrix>::value_type;
    using index_type = typename SparseTraits<DokMatrix>::index_type;
    using reference = typename SparseTraits<DokMatrix>::reference;
    using const_reference = typename SparseTraits<DokMatrix>::const_reference;
    using size_type = typename SparseTraits<DokMatrix>::size_type;
    using difference_type = typename SparseTraits<DokMatrix>::difference_type;
    using allocator_type = typename SparseTraits<DokMatrix>::allocator_type;
    using dictionary_type = typename SparseTraits<DokMatrix>::dictionary_type;

    [[using gnu: always_inline]]
    explicit DokMatrix(const allocator_type& alloc = {}) noexcept
        : m_dictionary(typename dictionary_type::allocator_type{alloc}) {}
    DokMatrix(const DokMatrix& other) noexcept = default;
    auto operator=(const DokMatrix& other) noexcept -> DokMatrix& = default;
    DokMatrix(DokMatrix&& other) noexcept = default;
    auto operator=(DokMatrix&& other) noexcept -> DokMatrix& = default;
    ~DokMatrix() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_dictionary.get_allocator();
    }

    [[using gnu: always_inline]]
    explicit DokMatrix(size_type nrows, size_type ncols, const allocator_type& alloc = {}) noexcept
        : m_nrows{nrows}, m_ncols{ncols}, m_dictionary(alloc) {}

    [[using gnu: always_inline]]
    DokMatrix(
        std::initializer_list<std::pair<const IndexPair<Index>, value_type>> initialize_list,
        const allocator_type& alloc = {}
    ) noexcept
        : m_dictionary{std::move(initialize_list), alloc} {
        for (const auto& [index_pair, coeff] : m_dictionary) {
            m_nrows = std::max(m_nrows, static_cast<size_type>(index_pair.row + 1));
            m_ncols = std::max(m_ncols, static_cast<size_type>(index_pair.col + 1));
        }
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
        return m_dictionary.size();
    }

    [[using gnu: always_inline]]
    auto resize(size_type nrows, size_type ncols) noexcept -> void {
#ifndef BOYLE_USE_BOOST_UNORDERED
        using std::erase_if;
#else
        using boost::unordered::erase_if;
#endif
        erase_if(m_dictionary, [nrows, ncols](const auto& item) constexpr noexcept -> bool {
            return item.first.row >= static_cast<index_type>(nrows) ||
                   item.first.col >= static_cast<index_type>(ncols);
        });
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto reserve(size_type capacity) noexcept -> void {
        m_dictionary.reserve(capacity);
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_dictionary.clear();
        return;
    }

    [[using gnu: always_inline, leaf]]
    auto compress() noexcept -> void {
#ifndef BOYLE_USE_BOOST_UNORDERED
        using std::erase_if;
#else
        using boost::unordered::erase_if;
#endif
        erase_if(m_dictionary, [](const auto& item) constexpr noexcept -> bool {
            return item.second == static_cast<value_type>(0.0);
        });
        return;
    }

    [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
    auto coeff(size_type row, size_type col) const noexcept -> value_type {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return static_cast<value_type>(0.0);
        }
        if (const auto search =
                m_dictionary.find({static_cast<index_type>(row), static_cast<index_type>(col)});
            search != m_dictionary.cend()) {
            return search->second;
        }
        return static_cast<value_type>(0.0);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(size_type row, size_type col, const_reference value) noexcept -> void {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return;
        }
        m_dictionary[{static_cast<index_type>(row), static_cast<index_type>(col)}] = value;
        return;
    }

    [[using gnu: flatten, leaf, hot]]
    auto eraseCoeff(size_type row, size_type col) noexcept -> void {
        if (row >= m_nrows || col >= m_ncols) [[unlikely]] {
            return;
        }
        if (auto search =
                m_dictionary.find({static_cast<index_type>(row), static_cast<index_type>(col)});
            search != m_dictionary.end()) {
            m_dictionary.erase(search);
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type row, size_type col) const noexcept -> value_type {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto dictionary() const noexcept -> const dictionary_type& {
        return m_dictionary;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_nrows;
        archive & m_ncols;
        archive & m_dictionary;
        return;
    }

    size_type m_nrows{0};
    size_type m_ncols{0};
    dictionary_type m_dictionary{};
};

} // namespace boyle::math
