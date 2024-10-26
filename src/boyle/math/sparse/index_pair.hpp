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

#include <concepts>

#include "boost/container_hash/hash.hpp"

namespace boyle::math {

template <std::integral Index = int>
struct IndexPair final {
    Index row;
    Index col;
};

template <std::integral Index = int>
struct IndexPairHash final {
    [[using gnu: const, always_inline, hot]]
    constexpr auto operator()(const IndexPair<Index>& index_pair) const noexcept -> std::size_t {
        std::size_t seed = 0;
        boost::hash_combine(seed, index_pair.row);
        boost::hash_combine(seed, index_pair.col);
        return seed;
    }
};

template <std::integral Index = int>
struct IndexPairEqual final {
    [[using gnu: const, always_inline, leaf, hot]]
    constexpr auto operator()(const IndexPair<Index>& lhs, const IndexPair<Index>& rhs)
        const noexcept -> bool {
        return lhs.row == rhs.row && lhs.col == rhs.col;
    }
};

template <std::integral Index = int>
struct IndexPairRowMajorCompare final {
    [[using gnu: const, always_inline, leaf, hot]]
    constexpr auto operator()(const IndexPair<Index>& lhs, const IndexPair<Index>& rhs)
        const noexcept -> bool {
        return lhs.row != rhs.row ? lhs.row < rhs.row : lhs.col < rhs.col;
    }
};

template <std::integral Index = int>
struct IndexPairColumnMajorCompare final {
    [[using gnu: const, always_inline, leaf, hot]]
    constexpr auto operator()(const IndexPair<Index>& lhs, const IndexPair<Index>& rhs)
        const noexcept -> bool {
        return lhs.col != rhs.col ? lhs.col < rhs.col : lhs.row < rhs.row;
    }
};

} // namespace boyle::math

namespace boost::serialization {

template <std::integral Index = int>
[[using gnu: always_inline]]
inline constexpr auto serialize(
    auto& archive, boyle::math::IndexPair<Index>& index_pair,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & index_pair.row;
    archive & index_pair.col;
    return;
}

} // namespace boost::serialization
