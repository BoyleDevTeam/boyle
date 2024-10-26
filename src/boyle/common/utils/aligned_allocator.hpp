/**
 * @file aligned_allocator.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-03-24
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <cstdlib>
#include <limits>
#include <memory>

namespace boyle::common {

template <typename T, std::size_t Alignment>
class AlignedAllocator : public std::allocator<T> {
  public:
    using value_type = T;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <class U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    AlignedAllocator() noexcept = default;
    AlignedAllocator(const AlignedAllocator& other) noexcept = default;
    auto operator=(const AlignedAllocator& other) noexcept -> AlignedAllocator& = default;
    ~AlignedAllocator() noexcept {}

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>& other) noexcept
        : std::allocator<U>(other) {}

    [[nodiscard]]
    static constexpr auto max_size() noexcept -> size_type {
        return std::numeric_limits<std::ptrdiff_t>::max() / sizeof(T);
    }

    auto allocate(size_type num) -> pointer {
        return static_cast<pointer>(std::aligned_alloc(Alignment, num * sizeof(T)));
    }

    auto deallocate(pointer p, size_type /*num*/) noexcept -> void {
        std::free(p);
        return;
    }
};

} // namespace boyle::common
