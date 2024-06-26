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
#include <vector>

#include "Eigen/Eigen"

#include "boyle/common/utils/macros.hpp"

namespace boyle::common {

template <typename T>
class ArrayView final {
  public:
    using value_type = T;
    using iterator = std::add_pointer_t<T>;
    using const_iterator = std::add_const_t<iterator>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    DISABLE_IMPLICIT_CONSTRUCTORS(ArrayView);

    ~ArrayView() noexcept {
        if (m_allocated) {
            delete[] m_data;
        }
    };

    // NOLINTNEXTLINE(modernize-avoid-c-arrays)
    explicit ArrayView(T data[], std::size_t size) noexcept : m_data{data}, m_size{size} {}

    explicit ArrayView(std::initializer_list<T>&& other) noexcept
        : m_data{new T[other.size()]}, m_size{other.size()}, m_allocated{true} {
        std::copy(other.begin(), other.end(), m_data);
    }

    template <std::contiguous_iterator ContiguousIt>
    ArrayView(ContiguousIt first, std::sentinel_for<ContiguousIt> auto last) noexcept
        requires std::same_as<
                     std::remove_const_t<typename std::iterator_traits<ContiguousIt>::value_type>,
                     T>
        : m_data{const_cast<T*>(first)}, m_size{static_cast<std::size_t>(last - first)} {}

    template <std::size_t Size>
    ArrayView(const std::array<T, Size>& other) noexcept
        : m_data{const_cast<T*>(other.data())}, m_size{other.size()} {}

    template <std::size_t Size>
    ArrayView(std::array<T, Size>&& other) noexcept
        : m_data{new T[other.size()]}, m_size{other.size()}, m_allocated(true) {
        std::copy(other.cbegin(), other.cend(), m_data);
    }

    template <std::size_t Size>
    operator std::array<T, Size>() const noexcept {
        std::array<T, Size> std_array;
        std::copy(m_data, m_data + m_size, std_array.begin());
        return std_array;
    }

    ArrayView(const std::vector<T>& other) noexcept
        : m_data{const_cast<T*>(other.data())}, m_size{other.size()} {}

    ArrayView(std::vector<T>&& other) noexcept
        : m_data{new T[other.size()]}, m_size{other.size()}, m_allocated(true) {
        std::copy(other.cbegin(), other.cend(), m_data);
    }

    operator std::vector<T>() const noexcept {
        std::vector<T> std_vector;
        std_vector.reserve(m_size);
        std_vector.assign(m_data, m_data + m_size);
        return std_vector;
    }

    template <int Size>
    ArrayView(const Eigen::Vector<T, Size>& other) noexcept
        : m_data{const_cast<T*>(other.data())}, m_size{static_cast<std::size_t>(other.size())} {}

    template <int Size>
    ArrayView(Eigen::Vector<T, Size>&& other) noexcept
        : m_data{new T[other.size()]}, m_size{static_cast<std::size_t>(other.size())},
          m_allocated(true) {
        std::copy(other.cbegin(), other.cend(), m_data);
    }

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
        return m_data == other.m_data && m_size == other.m_size;
    }

  private:
    T* m_data;
    std::size_t m_size;
    bool m_allocated{false};
};

} // namespace boyle::common
