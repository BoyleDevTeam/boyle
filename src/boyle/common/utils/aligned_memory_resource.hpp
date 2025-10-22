/**
 * @file aligned_memory_resource.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-06-25
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <cstddef>
#include <cstdlib>
#include <memory_resource>
#include <new>

namespace boyle::common::pmr {

class [[nodiscard]] AlignedMemoryResource final : public std::pmr::memory_resource {
  public:
    [[using gnu: always_inline]]
    explicit AlignedMemoryResource(std::size_t alignment = 32)
        : m_alignment{alignment} {
        if (alignment % 2 != 0) {
            throw std::invalid_argument("Alignment must be a power of 2.");
        }
    }
    constexpr AlignedMemoryResource(const AlignedMemoryResource& other) noexcept = default;
    constexpr auto operator=(const AlignedMemoryResource& other) noexcept
        -> AlignedMemoryResource& = default;
    constexpr AlignedMemoryResource(AlignedMemoryResource&& other) noexcept = default;
    constexpr auto operator=(AlignedMemoryResource&& other) noexcept
        -> AlignedMemoryResource& = default;
    ~AlignedMemoryResource() noexcept override = default;

  protected:
    [[using gnu: pure, always_inline, hot]]
    auto do_allocate(std::size_t bytes, [[maybe_unused]] std::size_t alignment) -> void* override {
        void* ptr{std::aligned_alloc(m_alignment, bytes)};
        if (ptr == nullptr) [[unlikely]] {
            throw std::bad_alloc();
        }
        return ptr;
    }

    [[using gnu: always_inline, hot]]
    auto do_deallocate(
        void* ptr, [[maybe_unused]] std::size_t bytes, [[maybe_unused]] std::size_t alignment
    ) -> void override {
        std::free(ptr);
        return;
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    auto do_is_equal(const std::pmr::memory_resource& other) const noexcept -> bool override {
        return this == &other;
    }

  private:
    std::size_t m_alignment;
};

[[using gnu: always_inline, hot]]
inline constexpr auto getAlignedMemoryResource(std::size_t alignment = 32) noexcept
    -> AlignedMemoryResource* {
    static AlignedMemoryResource memory_resource(alignment);
    return &memory_resource;
}

} // namespace boyle::common::pmr
