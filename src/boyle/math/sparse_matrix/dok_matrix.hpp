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
#include "boyle/math/sparse_matrix/index_pair.hpp"

namespace boyle::math {

template <GeneralArithmetic Scalar = double, std::integral Index = int>
class DokMatrix final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    using index_type = Index;

    DokMatrix() noexcept = default;
    DokMatrix(const DokMatrix& other) noexcept = default;
    auto operator=(const DokMatrix& other) noexcept -> DokMatrix& = default;
    DokMatrix(DokMatrix&& other) noexcept = default;
    auto operator=(DokMatrix&& other) noexcept -> DokMatrix& = default;
    ~DokMatrix() noexcept = default;

    [[using gnu: always_inline]]
    explicit DokMatrix(std::size_t nrows, std::size_t ncols) noexcept
        : m_nrows{nrows}, m_ncols{ncols} {}

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
        return m_dictionary.size();
    }

    [[using gnu: always_inline]]
    auto resize(std::size_t nrows, std::size_t ncols) noexcept -> void {
        boost::unordered::erase_if(m_dictionary, [nrows, ncols](const auto& item) noexcept -> bool {
            return item.first.row >= static_cast<Index>(nrows) ||
                   item.first.col >= static_cast<Index>(ncols);
        });
        m_nrows = nrows;
        m_ncols = ncols;
        return;
    }

    [[using gnu: always_inline]]
    auto reserve(std::size_t capacity) noexcept -> void {
        m_dictionary.reserve(capacity);
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_dictionary.clear();
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
        if (const auto search = m_dictionary.find({row, col}); search != m_dictionary.cend()) {
            return search->second;
        }
        return 0.0;
    }

    [[using gnu: flatten, leaf, hot]]
    auto updateCoeff(Index row, Index col, Scalar value) noexcept -> void {
        if (row >= static_cast<Index>(m_nrows) || col >= static_cast<Index>(m_ncols)) [[unlikely]] {
            return;
        }
        if (auto search = m_dictionary.find({row, col}); search != m_dictionary.end()) {
            if (value != 0.0) {
                search->second = value;
            } else {
                m_dictionary.erase(search);
            }
        } else {
            if (value != 0.0) {
                m_dictionary.emplace(IndexPair<Index>{row, col}, value);
            }
        }
        return;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(Index row, Index col) const noexcept -> Scalar {
        return coeff(row, col);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto dictionary() const noexcept
        -> const boost::unordered_flat_map<
            IndexPair<Index>, Scalar, IndexPairHash<Index>, IndexPairEqual<Index>>& {
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

    std::size_t m_nrows{0};
    std::size_t m_ncols{0};
    boost::unordered_flat_map<IndexPair<Index>, Scalar, IndexPairHash<Index>, IndexPairEqual<Index>>
        m_dictionary{};
};

} // namespace boyle::math
