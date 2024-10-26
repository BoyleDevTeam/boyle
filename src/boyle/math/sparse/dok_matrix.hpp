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
#include "boost/unordered/unordered_flat_map.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/sparse/index_pair.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar = double, std::integral Index = int>
class DokMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    DokMatrix() noexcept = default;
    DokMatrix(const DokMatrix& other) noexcept = default;
    auto operator=(const DokMatrix& other) noexcept -> DokMatrix& = default;
    DokMatrix(DokMatrix&& other) noexcept = default;
    auto operator=(DokMatrix&& other) noexcept -> DokMatrix& = default;
    ~DokMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit DokMatrix(size_type nrows, size_type ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

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
        boost::unordered::erase_if(m_dictionary, [nrows, ncols](const auto& item) noexcept -> bool {
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
    }

    [[using gnu: always_inline, leaf]]
    auto compress() noexcept -> void {
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
        if (auto search =
                m_dictionary.find({static_cast<index_type>(row), static_cast<index_type>(col)});
            search != m_dictionary.end()) {
            if (value != static_cast<value_type>(0.0)) {
                search->second = value;
            } else {
                m_dictionary.erase(search);
            }
        } else {
            if (value != static_cast<value_type>(0.0)) {
                m_dictionary.emplace(
                    IndexPair<index_type>{
                        static_cast<index_type>(row), static_cast<index_type>(col)
                    },
                    value
                );
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto operator[](size_type row, size_type col) const noexcept -> value_type {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto dictionary() const noexcept -> const boost::unordered_flat_map<
        IndexPair<index_type>, Scalar, IndexPairHash<index_type>, IndexPairEqual<index_type>>& {
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
    boost::unordered_flat_map<
        IndexPair<index_type>, value_type, IndexPairHash<index_type>, IndexPairEqual<index_type>>
        m_dictionary{};
};

} // namespace boyle::math
