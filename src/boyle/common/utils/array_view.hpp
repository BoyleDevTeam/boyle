/**
 * @file array_view.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-22
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <format>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Eigen/Eigen"

#include "boyle/common/utils/macros.hpp"

namespace boyle::common {

template <typename T, typename Allocator = std::allocator<T>>
class ArrayView final {
  public:
    using value_type = T;
    using allocator_type = Allocator;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;

    using iterator = std::add_pointer_t<T>;
    using const_iterator = std::add_const_t<iterator>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    DISABLE_IMPLICIT_CONSTRUCTORS(ArrayView);

    ~ArrayView() noexcept {
        if (m_allocated) {
            m_allocator.deallocate(m_data, m_size);
        }
    };

    explicit ArrayView(std::size_t size)
        : m_allocated{true}, m_size{size}, m_data{m_allocator.allocate(m_size)} {}

    // NOLINTNEXTLINE(modernize-avoid-c-arrays)
    explicit ArrayView(std::size_t size, T data[]) noexcept : m_size{size}, m_data{data} {}

    explicit ArrayView(std::initializer_list<T>&& other)
        : m_allocated{true}, m_size{other.size()}, m_data{m_allocator.allocate(m_size)} {
        std::copy(other.begin(), other.end(), m_data);
    }

    template <std::contiguous_iterator ContiguousIt>
    ArrayView(ContiguousIt first, std::sentinel_for<ContiguousIt> auto last) noexcept
        requires std::same_as<
                     std::remove_const_t<typename std::iterator_traits<ContiguousIt>::value_type>,
                     T>
        : m_size{static_cast<std::size_t>(last - first)}, m_data{const_cast<T*>(first)} {}

    template <std::size_t Size>
    ArrayView(const std::array<T, Size>& other) noexcept
        : m_size{other.size()}, m_data{const_cast<T*>(other.data())} {}

    template <std::size_t Size>
    ArrayView(std::array<T, Size>&& other) noexcept : m_size{other.size()}, m_data{other.data()} {}

    template <std::size_t Size>
    operator std::array<T, Size>() const noexcept {
        std::array<T, Size> std_array;
        std::copy(m_data, m_data + m_size, std_array.begin());
        return std_array;
    }

    ArrayView(const std::vector<T>& other) noexcept
        : m_size{other.size()}, m_data{const_cast<T*>(other.data())} {}

    ArrayView(std::vector<T>&& other) noexcept : m_size{other.size()}, m_data{other.data()} {}

    operator std::vector<T>() const noexcept {
        std::vector<T> std_vector;
        std_vector.reserve(m_size);
        std_vector.assign(m_data, m_data + m_size);
        return std_vector;
    }

    template <int Size>
    ArrayView(const Eigen::Vector<T, Size>& other) noexcept
        : m_size{static_cast<std::size_t>(other.size())}, m_data{const_cast<T*>(other.data())} {}

    template <int Size>
    ArrayView(Eigen::Vector<T, Size>&& other) noexcept
        : m_size{static_cast<std::size_t>(other.size())}, m_data{other.data()} {}

    template <int Size>
    operator Eigen::Vector<T, Size>() const noexcept {
        Eigen::Vector<T, Size> eigen_vector;
        eigen_vector.resize(m_size);
        std::copy(m_data, m_data + m_size, eigen_vector.data());
        return eigen_vector;
    }

    auto at(std::size_t pos) -> T& {
        if (pos >= m_size) {
            throw std::out_of_range(std::format(
                "Out-of-range error detected! The input pos has to be in the allowed range: pos = "
                "{0:d} while max_pos = {1:d}).",
                pos, m_size
            ));
        }
        return *(m_data + pos);
    }

    auto at(std::size_t pos) const -> const T& {
        if (pos >= m_size) {
            throw std::out_of_range(std::format(
                "Out-of-range error detected! The input pos has to be in the allowed range: pos = "
                "{0:d} while max_pos = {1:d}).",
                pos, m_size
            ));
        }
        return *(m_data + pos);
    }

    auto operator[](std::size_t pos) noexcept -> T& { return *(m_data + pos); }
    auto operator[](std::size_t pos) const noexcept -> const T& { return *(m_data + pos); }
    auto front() noexcept -> T& { return *m_data; }
    auto back() noexcept -> T& { return *(m_data + m_size - 1); }
    auto front() const noexcept -> const T& { return *m_data; }
    auto back() const noexcept -> const T& { return *(m_data + m_size - 1); }
    auto data() noexcept -> T* { return m_data; }
    auto data() const noexcept -> const T* { return m_data; }

    auto begin() noexcept -> iterator { return static_cast<iterator>(m_data); }
    auto end() noexcept -> iterator { return static_cast<iterator>(m_data + m_size); }
    auto cbegin() const noexcept -> const_iterator { return static_cast<const_iterator>(m_data); }
    auto cend() const noexcept -> const_iterator {
        return static_cast<const_iterator>(m_data + m_size);
    }
    auto rbegin() noexcept -> reverse_iterator {
        return std::make_reverse_iterator(m_data + m_size);
    }
    auto rback() noexcept -> reverse_iterator { return std::make_reverse_iterator(m_data); }
    auto crbegin() const noexcept -> const_reverse_iterator {
        return std::make_reverse_iterator<const_reverse_iterator>(m_data + m_size);
    }
    auto crback() const noexcept -> const_reverse_iterator {
        return std::make_reverse_iterator<const_reverse_iterator>(m_data);
    }

    [[nodiscard]]
    auto empty() const noexcept -> bool {
        return m_size == 0;
    }
    [[nodiscard]]
    auto size() const noexcept -> std::size_t {
        return m_size;
    }

    auto operator==(const ArrayView& other) const noexcept -> bool {
        return m_size == other.m_size && m_data == other.m_data;
    }

  private:
    [[no_unique_address]] allocator_type m_allocator{};

    bool m_allocated{false};
    std::size_t m_size{0};
    pointer m_data{nullptr};
};

} // namespace boyle::common

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

template <typename T>
struct is_array<::boyle::common::ArrayView<T>> final : public true_type {};

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
