/**
 * @file array.hpp
 * @author Houchen Li (houchenli@outlook.com)
 * @brief
 * @version 0.1
 * @date 2024-07-15
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved..
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <format>
#include <initializer_list>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

template <std::floating_point T, typename Allocator = std::allocator<T>>
class Array final {
    friend class boost::serialization::access;

  public:
    using element_type = T;
    using value_type = std::remove_cv_t<element_type>;
    using allocator_type = Allocator;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;

    using iterator = value_type*;
    using const_iterator = const value_type*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    [[using gnu: always_inline]]
    Array() noexcept = default;
    [[using gnu: always_inline]]
    Array(const Array& other) noexcept = default;
    [[using gnu: always_inline]]
    Array(Array&& other) noexcept = default;
    [[using gnu: always_inline]]
    auto operator=(const Array& other) noexcept -> Array& = default;
    [[using gnu: always_inline]]
    auto operator=(Array&& other) noexcept -> Array& = default;
    [[using gnu: always_inline]]
    ~Array() noexcept = default;

    [[using gnu: always_inline]]
    explicit Array(std::size_t size) noexcept
        : m_data(size) {
        m_data.shrink_to_fit();
    }
    [[using gnu: always_inline]]
    explicit Array(std::size_t size, const value_type& value) noexcept
        : m_data(size, value) {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    Array(std::initializer_list<value_type> obj) noexcept
        : m_data{std::move(obj)} {
        m_data.shrink_to_fit();
    }

    [[using gnu: always_inline]]
    Array(const InstanceOfTemplate<Array> auto& other) noexcept
        : m_data(other.size()) {
        std::transform(
            other.cbegin(), other.cend(), begin(),
            [](const decltype(other)::value_type& a) noexcept -> value_type {
                return static_cast<value_type>(a);
            }
        );
    }

    [[using gnu: always_inline]]
    Array(std::vector<value_type, allocator_type> obj) noexcept
        : m_data{std::move(obj)} {
        m_data.shrink_to_fit();
    }

    [[using gnu: pure,
      always_inline]] operator std::vector<value_type, allocator_type>() const& noexcept {
        return m_data;
    }
    [[using gnu: always_inline]] operator std::vector<value_type, allocator_type>() && noexcept {
        return std::move(m_data);
    }

    [[using gnu: always_inline]]
    auto assign(std::size_t size) noexcept -> void {
        m_data.assign(size);
        m_data.shrink_to_fit();
    }
    [[using gnu: always_inline]]
    auto assign(std::size_t size, const value_type& value) noexcept -> void {
        m_data.assign(size, value);
        m_data.shrink_to_fit();
    }

    [[using gnu: pure, always_inline]]
    auto at(std::size_t pos) -> reference {
        return m_data.at(pos);
    }
    [[using gnu: pure, always_inline]]
    auto at(std::size_t pos) const -> const_reference {
        return m_data.at(pos);
    }
    [[using gnu: pure, always_inline]]
    auto operator[](std::size_t pos) noexcept -> reference {
        return m_data[pos];
    }
    [[using gnu: pure, always_inline]]
    auto operator[](std::size_t pos) const noexcept -> const_reference {
        return m_data[pos];
    }

    [[using gnu: pure, always_inline]]
    auto front() noexcept -> reference {
        return m_data.front();
    }
    [[using gnu: pure, always_inline]]
    auto back() noexcept -> reference {
        return m_data.back();
    }
    [[using gnu: pure, always_inline]]
    auto front() const noexcept -> const_reference {
        return m_data.front();
    }
    [[using gnu: pure, always_inline]]
    auto back() const noexcept -> const_reference {
        return m_data.back();
    }
    [[using gnu: pure, always_inline]]
    auto data() noexcept -> pointer {
        return m_data.data();
    }
    [[using gnu: pure, always_inline]]
    auto data() const noexcept -> const_pointer {
        return m_data.data();
    }

    [[using gnu: pure, always_inline]]
    auto begin() noexcept -> iterator {
        return m_data.data();
    }
    [[using gnu: pure, always_inline]]
    auto end() noexcept -> iterator {
        return m_data.data() + m_data.size();
    }
    [[using gnu: pure, always_inline]]
    auto cbegin() const noexcept -> const_iterator {
        return m_data.data();
    }
    [[using gnu: pure, always_inline]]
    auto cend() const noexcept -> const_iterator {
        return m_data.data() + m_data.size();
    }
    [[using gnu: pure, always_inline]]
    auto rbegin() noexcept -> reverse_iterator {
        return std::make_reverse_iterator(end());
    }
    [[using gnu: pure, always_inline]]
    auto rend() noexcept -> reverse_iterator {
        return std::make_reverse_iterator(begin());
    }
    [[using gnu: pure, always_inline]]
    auto crbegin() const noexcept -> const_reverse_iterator {
        return std::make_reverse_iterator<const_reverse_iterator>(cend());
    }
    [[using gnu: pure, always_inline]]
    auto crend() const noexcept -> const_reverse_iterator {
        return std::make_reverse_iterator<const_reverse_iterator>(cbegin());
    }

    [[nodiscard]] [[using gnu: pure, always_inline]]
    auto empty() const noexcept -> bool {
        return m_data.empty();
    }
    [[nodiscard]] [[using gnu: pure, always_inline]]
    auto size() const noexcept -> std::size_t {
        return m_data.size();
    }

    [[using gnu: always_inline]]
    auto operator+=(const InstanceOfTemplate<Array> auto& obj) & noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array& {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::operator \"+\" requires the sizes "
                "of two arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        std::transform(cbegin(), cend(), obj.cbegin(), begin(), std::plus<value_type>{});
        return *this;
    }

    [[using gnu: pure, always_inline]]
    auto operator+(const InstanceOfTemplate<Array> auto& obj) const& noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array {
        Array result{*this};
        result += obj;
        return result;
    }

    [[using gnu: always_inline]]
    auto operator+(const InstanceOfTemplate<Array> auto& obj) && noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array&& {
        this->operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline]]
    auto operator-=(const InstanceOfTemplate<Array> auto& obj) & noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array& {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::operator \"-\" requires the sizes "
                "of two arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        std::transform(cbegin(), cend(), obj.cbegin(), begin(), std::minus<value_type>{});
        return *this;
    }

    [[using gnu: pure, always_inline]]
    auto operator-(const InstanceOfTemplate<Array> auto& obj) const& noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array {
        Array result{*this};
        result -= obj;
        return result;
    }

    [[using gnu: always_inline]]
    auto operator-(const InstanceOfTemplate<Array> auto& obj) && noexcept(!BOYLE_CHECK_PARAMS
    ) -> Array&& {
        this->operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline]]
    auto operator*=(Arithmetic auto fac) & noexcept(!BOYLE_CHECK_PARAMS) -> Array& {
#if BOYLE_CHECK_PARAMS == 1
        if (std::isnan(fac)) {
            throw std::invalid_argument(
                "Invalid arguments detected! boyle::math::Array::operator \"*=\" does not accept "
                "NAN as a factor."
            );
        }
#endif
        std::transform(
            cbegin(), cend(), begin(),
            [fac](const value_type& x) noexcept -> value_type { return x * fac; }
        );
        return *this;
    }

    [[using gnu: pure, always_inline]]
    auto operator*(Arithmetic auto fac) const& noexcept(!BOYLE_CHECK_PARAMS) -> Array {
        Array result{*this};
        result *= fac;
        return result;
    }

    [[using gnu: always_inline]]
    auto operator*(Arithmetic auto fac) && noexcept(!BOYLE_CHECK_PARAMS) -> Array&& {
        this->operator*=(fac);
        return std::move(*this);
    }

    [[using gnu: always_inline]]
    auto operator/=(Arithmetic auto den) & noexcept(!BOYLE_CHECK_PARAMS) -> Array& {
#if BOYLE_CHECK_PARAMS == 1
        if (den == 0.0) {
            throw std::invalid_argument(
                "Invalid arguments detected! boyle::math::Array::operator \"/=\" does not accept "
                "0.0 as a denominator."
            );
        }
#endif
        const value_type fac = 1.0 / den;
        std::transform(
            cbegin(), cend(), begin(),
            [fac](const value_type& x) noexcept -> value_type { return x * fac; }
        );
        return *this;
    }

    [[using gnu: pure, always_inline]]
    auto operator/(Arithmetic auto den) const& noexcept(!BOYLE_CHECK_PARAMS) -> Array {
        Array result{*this};
        result /= den;
        return result;
    }

    [[using gnu: always_inline]]
    auto operator/(Arithmetic auto den) && noexcept(!BOYLE_CHECK_PARAMS) -> Array&& {
        this->operator/=(den);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline]]
    auto operator-() const noexcept(!BOYLE_CHECK_PARAMS) -> Array {
        return this->operator*(-1.0);
    }

    [[using gnu: pure, always_inline]]
    auto dot(const InstanceOfTemplate<Array> auto& obj) const
        noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::dot() requires the sizes of two "
                "arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        return std::transform_reduce(cbegin(), cend(), obj.cbegin(), 0.0);
    }

    [[using gnu: pure, always_inline]]
    auto norm() const noexcept -> value_type {
        return std::sqrt(std::transform_reduce(cbegin(), cend(), cbegin(), 0.0));
    }

    [[using gnu: pure, always_inline]]
    auto euclideanTo(const InstanceOfTemplate<Array> auto& obj) const
        noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::euclideanTo() requires the sizes "
                "of two arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        return this->operator-(obj).norm();
    }

    [[using gnu: pure, always_inline]]
    auto identicalTo(const InstanceOfTemplate<Array> auto& obj, std::floating_point auto tol = 1E-8)
        const noexcept(!BOYLE_CHECK_PARAMS) -> bool {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::identicalTo() requires the sizes "
                "of two arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        return std::transform_reduce(
                   cbegin(), cend(), obj.cbegin(), 0.0,
                   [](const value_type& a, const value_type& b) noexcept -> value_type {
                       return std::max(a, b);
                   },
                   [](const value_type& a, const value_type& b) noexcept -> value_type {
                       return std::abs(a - b) / std::max(std::abs(a), 1.0);
                   }
               ) < tol;
    }

    [[using gnu: pure, always_inline]]
    auto orthogonalTo(
        const InstanceOfTemplate<Array> auto& obj, std::floating_point auto tol = 1E-8
    ) const noexcept(!BOYLE_CHECK_PARAMS) -> bool {
#if BOYLE_CHECK_PARAMS == 1
        if (size() != obj.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! boyle::math::Array::orthogonalTo() requires the sizes "
                "of two arrays to be identical: this->size() == {0:d} while obj.size() == {1:d}.",
                size(), obj.size()
            ));
        }
#endif
        const value_type test = this->dot(obj);
        return (test + tol) * (test - tol) < tol;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_data;
        return;
    }

    std::vector<value_type, allocator_type> m_data{};
};

[[using gnu: pure, always_inline]]
inline auto operator*(Arithmetic auto fac, const InstanceOfTemplate<Array> auto& obj) noexcept
    -> decltype(obj) {
    return obj * fac;
}

[[using gnu: always_inline]]
inline auto operator*(Arithmetic auto fac, InstanceOfTemplate<Array> auto&& obj) noexcept
    -> decltype(obj) {
    return std::move(obj) * fac;
}

[[using gnu: pure, always_inline]]
inline auto dot(
    const InstanceOfTemplate<Array> auto& obj_0, const InstanceOfTemplate<Array> auto& obj_1
) noexcept -> typename decltype(obj_0)::value_type {
    return obj_0.dot(obj_1);
}

[[using gnu: pure, always_inline]]
inline auto norm(const InstanceOfTemplate<Array> auto& obj) noexcept ->
    typename decltype(obj)::value_type {
    return obj.norm();
}

[[using gnu: pure, always_inline]]
inline auto euclidean(
    const InstanceOfTemplate<Array> auto& obj_0, const InstanceOfTemplate<Array> auto& obj_1
) noexcept -> typename decltype(obj_0)::value_type {
    return obj_0.euclideanTo(obj_1);
}

[[using gnu: pure, always_inline]]
inline auto identical(
    const InstanceOfTemplate<Array> auto& obj_0, const InstanceOfTemplate<Array> auto& obj_1,
    std::floating_point auto tol = 1E-8
) noexcept -> bool {
    return obj_0.identicalTo(obj_1, tol);
}

[[using gnu: pure, always_inline]]
inline auto orthogonal(
    const InstanceOfTemplate<Array> auto& obj_0, const InstanceOfTemplate<Array> auto& obj_1,
    std::floating_point auto tol = 1E-8
) noexcept -> bool {
    return obj_0.orthogonalTo(obj_1, tol);
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

inline auto hypot(const boyle::math::InstanceOfTemplate<boyle::math::Array> auto& obj) noexcept ->
    typename decltype(obj)::value_type {
    return obj.norm();
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
