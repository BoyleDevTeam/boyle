/**
 * @file dok_vector.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-06-19
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <initializer_list>

#include "boost/serialization/access.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/sparse_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar, std::integral Index, typename Alloc>
class DokVector final {
    friend class boost::serialization::access;

  public:
    using value_type = typename SparseTraits<DokVector>::value_type;
    using index_type = typename SparseTraits<DokVector>::index_type;
    using reference = typename SparseTraits<DokVector>::reference;
    using const_reference = typename SparseTraits<DokVector>::const_reference;
    using size_type = typename SparseTraits<DokVector>::size_type;
    using difference_type = typename SparseTraits<DokVector>::difference_type;
    using allocator_type = typename SparseTraits<DokVector>::allocator_type;
    using dictionary_type = typename SparseTraits<DokVector>::dictionary_type;

    [[using gnu: always_inline]]
    explicit DokVector(const allocator_type& alloc = {}) noexcept
        : m_dictionary{alloc} {}
    DokVector(const DokVector& other) noexcept = default;
    auto operator=(const DokVector& other) noexcept -> DokVector& = default;
    DokVector(DokVector&& other) noexcept = default;
    auto operator=(DokVector&& other) noexcept -> DokVector& = default;
    ~DokVector() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_dictionary.get_allocator();
    }

    [[using gnu: always_inline]]
    explicit DokVector(size_type n, const allocator_type& alloc = {}) noexcept
        : m_size(n, alloc) {}

    [[using gnu: always_inline]]
    DokVector(
        std::initializer_list<std::pair<const index_type, value_type>> initialize_list,
        const allocator_type& alloc = {}
    ) noexcept
        : m_dictionary{std::move(initialize_list), alloc} {
        for (const auto& [i, coeff] : m_dictionary) {
            m_size = std::max(m_size, static_cast<size_type>(i + 1));
        }
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto size() const noexcept -> size_type {
        return m_size;
    }

    [[using gnu: always_inline]]
    auto resize(size_type n) noexcept -> void {
#ifndef BOYLE_USE_BOOST_UNORDERED
        using std::erase_if;
#else
        using boost::unordered::erase_if;
#endif
        erase_if(m_dictionary, [n](const auto& item) constexpr noexcept -> bool {
            return static_cast<size_type>(item.first) >= n;
        });
        m_size = n;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto nnzs() const noexcept -> size_type {
        return m_dictionary.size();
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
    auto coeff(size_type i) const noexcept -> value_type {
        if (i >= m_size) [[unlikely]] {
            return static_cast<value_type>(0.0);
        }
        if (const auto search = m_dictionary.find(static_cast<index_type>(i));
            search != m_dictionary.cend()) {
            return search->second;
        }
        return static_cast<value_type>(0.0);
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(size_type i, const_reference value) noexcept -> void {
        if (i >= m_size) [[unlikely]] {
            return;
        }
        m_dictionary[static_cast<index_type>(i)] = value;
        return;
    }

    [[using gnu: flatten, leaf, hot]]
    auto eraseCoeff(size_type i) noexcept -> void {
        if (i >= m_size) [[unlikely]] {
            return;
        }
        if (auto search = m_dictionary.find(static_cast<index_type>(i));
            search != m_dictionary.end()) {
            m_dictionary.erase(search);
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type i) const noexcept -> value_type {
        return coeff(i);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto dictionary() const noexcept -> const dictionary_type& {
        return m_dictionary;
    }

  private:
    [[using gnu: always_inline]]
    void serialize(auto& ar, [[maybe_unused]] const unsigned int version) {
        ar & m_size;
        ar & m_dictionary;
        return;
    }

    size_type m_size{0};
    dictionary_type m_dictionary{};
};

} // namespace boyle::math
