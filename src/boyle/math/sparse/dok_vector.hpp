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

#ifndef BOYLE_USE_BOOST_UNORDERED
#include "boost/serialization/unordered_map.hpp"
#include <unordered_map>
#else
#include "boost/unordered/unordered_flat_map.hpp"
#endif

namespace boyle::math {

template <ScalarArithmetic Scalar = double, std::integral Index = int>
class DokVector final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using dictionary_type =
#ifndef BOYLE_USE_BOOST_UNORDERED
        std::unordered_map<index_type, value_type>;
#else
        boost::unordered::unordered_flat_map<index_type, value_type>;
#endif

    DokVector() noexcept = default;
    DokVector(const DokVector& other) noexcept = default;
    auto operator=(const DokVector& other) noexcept -> DokVector& = default;
    DokVector(DokVector&& other) noexcept = default;
    auto operator=(DokVector&& other) noexcept -> DokVector& = default;
    ~DokVector() noexcept = default;

    [[using gnu: always_inline]]
    explicit DokVector(size_type n) noexcept
        : m_size(n) {}

    [[using gnu: always_inline]]
    DokVector(std::initializer_list<std::pair<const index_type, value_type>>&& initialize_list
    ) noexcept
        : m_dictionary{std::move(initialize_list)} {
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
