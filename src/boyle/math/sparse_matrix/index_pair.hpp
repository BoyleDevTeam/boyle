/**
 * @file index_pair.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-19
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <bits/functional_hash.h>
#include <concepts>

namespace boyle::math {

template <std::integral Index = int>
struct IndexPair final {
    Index row;
    Index col;

    [[using gnu: pure, always_inline]]
    constexpr auto operator==(IndexPair index_pair) const noexcept -> bool {
        return row == index_pair.row && col == index_pair.col;
    }

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & row;
        archive & col;
        return;
    }
};

template <std::integral Index = int>
struct IndexPairHash final {
    [[using gnu: pure, always_inline]]
    constexpr auto operator()(IndexPair<Index> index_pair) const noexcept -> std::size_t {
        return std::hash<Index>()(index_pair.row) ^ std::hash<Index>()(index_pair.col);
    }
};

template <std::integral Index = int>
struct IndexPairRowMajorCompare final {
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto operator()(const IndexPair<Index>& lhs, const IndexPair<Index>& rhs)
        const noexcept -> bool {
        if (lhs.row < rhs.row) {
            return true;
        }
        if (lhs.row > rhs.row) {
            return false;
        }
        return lhs.col < rhs.col;
    }
};

template <std::integral Index = int>
struct IndexPairColumnMajorCompare final {
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto operator()(const IndexPair<Index>& lhs, const IndexPair<Index>& rhs)
        const noexcept -> bool {
        if (lhs.col < rhs.col) {
            return true;
        }
        if (lhs.col > rhs.col) {
            return false;
        }
        return lhs.row < rhs.row;
    }
};

} // namespace boyle::math
