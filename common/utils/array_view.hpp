/**
 * @file array_view.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-22
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "Eigen/Core"
#include "spdlog/spdlog.h"

#include "common/utils/tags.hpp"

namespace tiny_pnc {
namespace common {

template <typename T>
class ArrayView final {
  public:
    using value_type = T;
    using iterator = std::add_pointer_t<T>;
    using const_iterator = std::add_const_t<iterator>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    ArrayView() noexcept = default;
    ArrayView(const ArrayView& other) noexcept
        : data_(other.data_), size_(other.size_), allocated_(false) {}
    ArrayView(ArrayView&& other) noexcept
        : data_(other.data_), size_(other.size_), allocated_(false) {}
    ArrayView& operator=(const ArrayView& other) noexcept {
        data_ = other.data_;
        size_ = other.size_;
        allocated_ = false;
        return *this;
    }
    ArrayView& operator=(ArrayView&& other) noexcept {
        data_ = other.data_;
        size_ = other.size_;
        allocated_ = false;
        return *this;
    }
    ~ArrayView() noexcept {
        if (allocated_) {
            delete[] data_;
            data_ = nullptr;
            size_ = 0;
            allocated_ = false;
        }
    };

    explicit ArrayView(std::size_t size) noexcept : ArrayView(size, allocate_tag{}) {}

    explicit ArrayView(std::size_t size, allocate_tag tag) noexcept
        : data_(size ? new T[size] : nullptr), size_(size), allocated_(size ? true : false) {
        static_assert(!std::is_const_v<T>, "Should not allocate empty space for a const pointer.");
    }

    explicit ArrayView(T* data, std::size_t size) noexcept
        : data_(size ? data : nullptr), size_(size), allocated_(false) {}

    explicit ArrayView(T* data, std::size_t size, allocate_tag tag) noexcept {
        if (size) {
            size_ = size;
            std::remove_const_t<T>* temp_data = new std::remove_const_t<T>[size_];
            std::copy(data, data + size, temp_data);
            data_ = const_cast<T*>(temp_data);
            allocated_ = true;
        }
    }

    template <typename ForwardIt>
    ArrayView(ForwardIt begin, ForwardIt end) noexcept {
        static_assert(
            std::is_same_v<
                typename std::iterator_traits<ForwardIt>::value_type, std::remove_const_t<T>>,
            "The ForwardIt::value_type should be the same as ArrayView::value_type."
        );
        if (begin < end) {
            data_ = const_cast<T*>(begin);
            size_ = end - begin;
        }
    }

    template <typename ForwardIt>
    explicit ArrayView(ForwardIt begin, ForwardIt end, allocate_tag tag) noexcept {
        static_assert(
            std::is_same_v<
                typename std::iterator_traits<ForwardIt>::value_type, std::remove_const_t<T>>,
            "The ForwardIt::value_type should be the same as ArrayView::value_type."
        );
        if (begin < end) {
            size_ = end - begin;
            std::remove_const_t<T>* temp_data = new std::remove_const_t<T>[size_];
            std::copy(begin, end, temp_data);
            data_ = const_cast<T*>(temp_data);
            allocated_ = true;
        }
    }

    template <std::size_t Size>
    ArrayView(const std::array<std::remove_const_t<T>, Size>& other) noexcept {
        if (Size != 0) {
            data_ = const_cast<T*>(other.data());
            size_ = other.size();
        }
    }

    template <std::size_t Size>
    explicit ArrayView(
        const std::array<std::remove_const_t<T>, Size>& other, allocate_tag tag
    ) noexcept {
        if (other.size() != 0) {
            size_ = other.size();
            std::remove_const_t<T>* temp_data = new std::remove_const_t<T>[size_];
            std::copy(other.cbegin(), other.cend(), temp_data);
            data_ = const_cast<T*>(temp_data);
            allocated_ = true;
        }
    }

    template <std::size_t Size>
    operator std::array<std::remove_const_t<T>, Size>() {
        if (Size != size_) {
            std::string error_msg = std::format(
                "Out-of-range error detected! The initialized std::array size must equal to "
                "array_view.size(): std::array.size() = {0:d} while array_view.size = {1:d}).",
                Size, size_
            );
            throw std::out_of_range(std::move(error_msg));
        }
        std::array<std::remove_const_t<T>, Size> arr;
        std::copy(cbegin(), cend(), arr.data());
        return arr;
    }

    ArrayView(const std::vector<std::remove_const_t<T>>& other) noexcept {
        if (other.size() != 0) {
            data_ = const_cast<T*>(other.data());
            size_ = other.size();
        }
    }

    explicit ArrayView(
        const std::vector<std::remove_const_t<T>>& other, allocate_tag tag
    ) noexcept {
        if (other.size() != 0) {
            size_ = other.size();
            std::remove_const_t<T>* temp_data = new std::remove_const_t<T>[size_];
            std::copy(other.cbegin(), other.cend(), temp_data);
            data_ = const_cast<T*>(temp_data);
            allocated_ = true;
        }
    }

    operator std::vector<std::remove_const_t<T>>() noexcept {
        return std::vector<std::remove_const_t<T>>(cbegin(), cend());
    }

    template <int Size>
    ArrayView(const Eigen::Vector<std::remove_const_t<T>, Size>& other) noexcept {
        if (other.size()) {
            data_ = const_cast<T*>(other.data());
            size_ = other.size();
        }
    }

    template <int Size>
    explicit ArrayView(
        const Eigen::Vector<std::remove_const_t<T>, Size>& other, allocate_tag tag
    ) noexcept {
        if (other.size() != 0) {
            size_ = other.size();
            std::remove_const_t<T>* temp_data = new std::remove_const_t<T>[size_];
            std::copy(other.cbegin(), other.cend(), temp_data);
            data_ = const_cast<T*>(temp_data);
            allocated_ = true;
        }
    }

    template <int Size>
    operator Eigen::Vector<std::remove_const_t<T>, Size>() {
        if (Size != static_cast<int>(size_) && Size != Eigen::Dynamic) {
            std::string error_msg = std::format(
                "Out-of-range error detected! The initialized Eigen::Vector size must equal to "
                "ArrayView::size(): Eigen::Vector::size() = {0:d} while array_view.size = {1:d}).",
                Size, size_
            );
            throw std::out_of_range(std::move(error_msg));
        }
        Eigen::Vector<std::remove_const_t<T>, Size> arr(size_);
        std::copy(cbegin(), cend(), arr.data());
        return arr;
    }

    T& at(std::size_t pos) {
        if (pos >= size_) {
            std::string error_msg = std::format(
                "Out-of-range error detected! The input pos has to be in the allowed range: pos = "
                "{0:d} while max_pos = {1:d}).",
                pos, size_
            );
            throw std::out_of_range(std::move(error_msg));
        }
        return data_[pos];
    }

    const T& at(std::size_t pos) const {
        if (pos >= size_) {
            std::string error_msg = std::format(
                "Out-of-range error detected! The input pos has to be in the allowed range: pos = "
                "{0:d} while max_pos = {1:d}).",
                pos, size_
            );
            throw std::out_of_range(std::move(error_msg));
        }
        return data_[pos];
    }

    T& operator[](std::size_t pos) noexcept { return data_[pos]; }
    const T& operator[](std::size_t pos) const noexcept { return data_[pos]; }
    T& front() noexcept { return data_[0]; }
    T& back() noexcept { return data_[size_ - 1]; }
    const T& front() const noexcept { return data_[0]; }
    const T& back() const noexcept { return data_[size_ - 1]; }
    T* data() noexcept { return data_; }
    const T* data() const noexcept { return const_cast<const T*>(data); }

    iterator begin() noexcept { return static_cast<iterator>(data_); }
    iterator end() noexcept { return static_cast<iterator>(data_ + size_); }
    const_iterator cbegin() const noexcept { return static_cast<const_iterator>(data_); }
    const_iterator cend() const noexcept { return static_cast<const_iterator>(data_ + size_); }
    reverse_iterator rbegin() noexcept { return std::make_reverse_iterator(data_ + size_); }
    reverse_iterator rback() noexcept { return std::make_reverse_iterator(data_); }
    const_reverse_iterator crbegin() const noexcept {
        return std::make_reverse_iterator<const_reverse_iterator>(data_ + size_);
    }
    const_reverse_iterator crback() const noexcept {
        return std::make_reverse_iterator<const_reverse_iterator>(data_);
    }

    bool empty() const noexcept { return size_ == 0; }
    std::size_t size() const noexcept { return size_; }

    bool operator==(const ArrayView& other) const noexcept {
        if (data_ != other.data_ || size_ != other.size_ || allocated_ != other.allocated_) {
            return false;
        }
        return true;
    }

  private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
    bool allocated_ = false;
};

} // namespace common
} // namespace tiny_pnc
