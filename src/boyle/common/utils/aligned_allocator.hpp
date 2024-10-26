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

#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <type_traits>

namespace boyle::common {

template <typename T, std::size_t Alignment>
class AlignedAllocator {
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

    constexpr AlignedAllocator() noexcept = default;
    constexpr AlignedAllocator(const AlignedAllocator& other) noexcept = default;
    constexpr auto operator=(const AlignedAllocator& other) noexcept -> AlignedAllocator& = default;
    constexpr AlignedAllocator(AlignedAllocator&& other) noexcept = default;
    constexpr auto operator=(AlignedAllocator&& other) noexcept -> AlignedAllocator& = default;
    virtual ~AlignedAllocator() noexcept = default;

    template <typename U>
    constexpr AlignedAllocator([[maybe_unused]] const AlignedAllocator<U, Alignment>& other
    ) noexcept {}

    template <typename U, std::size_t OtherAlignment>
    constexpr auto operator==([[maybe_unused]] const AlignedAllocator<U, OtherAlignment>& other
    ) const noexcept -> bool {
        return std::is_same_v<T, U> && Alignment == OtherAlignment;
    }

    [[using gnu: always_inline]] [[nodiscard]]
    static constexpr auto max_size() noexcept -> size_type {
        return std::numeric_limits<std::ptrdiff_t>::max() / sizeof(T);
    }

    [[using gnu: pure, always_inline, hot]]
    auto allocate(size_type num) -> pointer {
        void* ptr{std::aligned_alloc(Alignment, num * sizeof(T))};
        if (ptr == nullptr) [[unlikely]] {
            throw std::bad_alloc();
        }
        return static_cast<pointer>(ptr);
    }

    [[using gnu: always_inline, hot]]
    auto deallocate(pointer p, size_type /*num*/) noexcept -> void {
        std::free(p);
        return;
    }

    auto select_on_container_copy_construction() const noexcept -> AlignedAllocator {
        return *this;
    }
};

} // namespace boyle::common
