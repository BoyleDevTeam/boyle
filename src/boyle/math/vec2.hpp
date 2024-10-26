/**
 * @file vec2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-17
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <cmath>
#include <concepts>
#include <format>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <vector>

#include "fmt/format.h"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

template <std::floating_point T>
struct Vec2 final {
    using value_type = T;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    [[using gnu: always_inline]]
    Vec2() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2(value_type cv) noexcept
        : x{cv}, y{cv} {}
    [[using gnu: always_inline]]
    constexpr Vec2(value_type cx, value_type cy) noexcept
        : x{cx}, y{cy} {}
    [[using gnu: always_inline]]
    constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Vec2& other) noexcept -> Vec2& = default;
    [[using gnu: always_inline]]
    constexpr Vec2(Vec2&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Vec2&& other) noexcept -> Vec2& = default;
    ~Vec2() noexcept = default;

    template <std::floating_point U>
    [[using gnu: always_inline]]
    constexpr explicit Vec2(const std::pair<U, U>& other) noexcept
        : x{static_cast<value_type>(other.first)}, y{static_cast<value_type>(other.second)} {}

    template <std::floating_point U>
    [[using gnu: always_inline]]
    constexpr explicit Vec2(const Vec2<U>& other) noexcept
        : x{static_cast<value_type>(other.x)}, y{static_cast<value_type>(other.y)} {}

    [[using gnu: always_inline]]
    constexpr auto operator[](size_t i) noexcept -> reference {
        return &x[i];
    }
    [[using gnu: always_inline]]
    constexpr auto operator[](size_t i) const noexcept -> const_reference {
        return &x[i];
    }

    [[using gnu: always_inline, leaf]]
    constexpr auto operator+=(const Vec2& obj) noexcept -> Vec2& {
        x += obj.x;
        y += obj.y;
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator+(const Vec2& obj) const& noexcept -> Vec2 {
        Vec2 result{*this};
        result += obj;
        return result;
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator+(const Vec2& obj) && noexcept -> Vec2&& {
        operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf]]
    constexpr auto operator-=(const Vec2& obj) noexcept -> Vec2& {
        x -= obj.x;
        y -= obj.y;
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator-(const Vec2& obj) const& noexcept -> Vec2 {
        Vec2 result{*this};
        result -= obj;
        return result;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator-(const Vec2& obj) && noexcept -> Vec2&& {
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf]]
    constexpr auto operator*=(const Arithmetic auto& fac) noexcept -> Vec2& {
        x *= static_cast<value_type>(fac);
        y *= static_cast<value_type>(fac);
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator*(const Arithmetic auto& fac) const& noexcept -> Vec2 {
        Vec2 result{*this};
        result *= fac;
        return result;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator*(const Arithmetic auto& fac) && noexcept -> Vec2&& {
        operator*=(fac);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf]]
    constexpr auto operator/=(const Arithmetic auto& den) noexcept -> Vec2& {
        operator*=(static_cast<decltype(den)>(1.0) / den);
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator/(const Arithmetic auto& den) const& noexcept -> Vec2 {
        Vec2 result{*this};
        result /= den;
        return result;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator/(const Arithmetic auto& den) && noexcept -> Vec2&& {
        operator/=(den);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-() const& noexcept -> Vec2 {
        return Vec2{-x, -y};
    }
    [[using gnu: pure, always_inline]]
    constexpr auto operator-() && noexcept -> Vec2&& {
        x = -x;
        y = -y;
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto operator==(const Vec2& other) const noexcept -> bool {
        return x == other.x && y == other.y;
    }

    template <std::floating_point U>
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(const Vec2<U>& obj) const noexcept -> value_type {
        return x * obj.x + y * obj.y;
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto crossProj(const Vec2<U>& obj) const noexcept -> value_type {
        return x * obj.y - y * obj.x;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto euclidean() const noexcept -> value_type {
        return std::hypot(x, y);
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanSqr() const noexcept -> value_type {
        return x * x + y * y;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto normalized() const noexcept -> Vec2 {
        return *this / euclidean();
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanTo(const Vec2<U>& obj) const noexcept -> value_type {
        return std::hypot(x - obj.x, y - obj.y);
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanSqrTo(const Vec2<U>& obj) const noexcept -> value_type {
        return (x - obj.x) * (x - obj.x) + (y - obj.y) * (y - obj.y);
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr auto identicalTo(const Vec2<U>& obj, value_type tol = 1E-8) const noexcept -> bool {
        return this->euclideanSqrTo(obj) < tol * tol;
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr auto orthogonalTo(const Vec2<U>& obj, value_type tol = 1E-8) const noexcept -> bool {
        return std::abs(this->dot(obj)) < tol;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto angle() const noexcept -> value_type {
        return std::atan2(y, x);
    }

    [[using gnu: pure, always_inline]]
    constexpr auto selfRotate(const std::floating_point auto& radian) noexcept -> Vec2& {
        *this = Vec2{
            x * std::cos(radian) - y * std::sin(radian), x * std::sin(radian) + y * std::cos(radian)
        };
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto rotate(const std::floating_point auto& radian) const& noexcept -> Vec2 {
        return Vec2{
            x * std::cos(radian) - y * std::sin(radian), x * std::sin(radian) + y * std::cos(radian)
        };
    }
    [[using gnu: pure, always_inline]]
    constexpr auto rotate(const std::floating_point auto& radian) && noexcept -> Vec2&& {
        this->selfRotate(radian);
        return std::move(*this);
    }

    [[using gnu: always_inline]]
    constexpr auto selfRotateHalfPi() noexcept -> Vec2& {
        *this = Vec2{-y, x};
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto rotateHalfPi() const& noexcept -> Vec2 {
        return Vec2{-y, x};
    }
    [[using gnu: pure, always_inline]]
    constexpr auto rotateHalfPi() && noexcept -> Vec2&& {
        this->selfRotateHalfPi();
        return std::move(*this);
    }

    value_type x;
    value_type y;
};

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;

template <typename Char, std::floating_point T>
[[using gnu: always_inline]]
inline auto operator<<(std::basic_ostream<Char>& os, const Vec2<T>& obj) noexcept
    -> std::basic_ostream<Char>& {
    os << "(x: " << obj.x << ", y: " << obj.y << ")";
    return os;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto operator*(const Arithmetic auto& fac, const Vec2<T>& obj) noexcept
    -> Vec2<T> {
    return obj * fac;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto operator*(const Arithmetic auto& fac, Vec2<T>&& obj) noexcept -> Vec2<T>&& {
    return std::move(obj) * fac;
}

template <InstanceOfTemplate<Vec2> T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline auto squeeze(
    std::ranges::forward_range auto&& xs, std::ranges::forward_range auto&& ys
) noexcept -> std::vector<T>
    requires std::same_as<std::ranges::range_value_t<decltype(xs)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(ys)>, typename T::value_type>
{
    if (xs.size() != ys.size()) {
        return std::vector<T>{};
    }
    std::vector<T> result;
    result.reserve(xs.size());
    for (auto it_xs{xs.begin()}, it_ys{ys.begin()}; !(it_xs == xs.end() || it_ys == ys.end());
         ++it_xs, ++it_ys) {
        result.emplace_back(*it_xs, *it_ys);
    }
    return result;
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

template <std::floating_point T, typename Char>
struct formatter<::boyle::math::Vec2<T>, Char> {
    [[using gnu: ]]
    constexpr auto parse(std::basic_format_parse_context<Char>& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin();
        auto end = ctx.end();
        if (it != end && *it >= '0' && *it <= '9') {
            width = 0;
            for (; it != end && *it >= '0' && *it <= '9'; ++it) {
                width = width * 10 + (*it - '0');
            }
        }
        if (it != end && *it == '.') {
            ++it;
            precision = 0;
            for (; it != end && *it >= '0' && *it <= '9'; ++it) {
                precision = precision * 10 + (*it - '0');
            }
        }
        return it;
    }

    template <typename FormatContext>
    [[using gnu: ]]
    auto format(const ::boyle::math::Vec2<T>& obj, FormatContext& ctx) const
        -> decltype(ctx.out()) {
        const std::uint8_t condition = ((width > 0) << 1) + (precision >= 0);
        std::basic_ostringstream<Char> ss;
        switch (condition) {
        case 1:
            ss << std::fixed << std::setprecision(precision) << "(x: " << obj.x << ", y: " << obj.y
               << ")";
            break;
        case 2:
            ss << std::fixed << "(x: " << std::setw(width) << obj.x << ", y: " << std::setw(width)
               << obj.y << ")";
            break;
        case 3:
            ss << std::fixed << std::setprecision(precision) << "(x: " << std::setw(width) << obj.x
               << ", y: " << std::setw(width) << obj.y << ")";
            break;
        default:
            ss << "(x: " << obj.x << ", y: " << obj.y << ")";
            break;
        }
        return std::format_to(ctx.out(), "{0:s}", ss.str());
    }

    int width{0};
    int precision{-1};
};

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto hypot(const ::boyle::math::Vec2<T>& obj) noexcept -> T {
    return std::hypot(obj.x, obj.y);
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto abs(const ::boyle::math::Vec2<T>& obj) noexcept -> T {
    return obj.euclidean();
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto norm(const ::boyle::math::Vec2<T>& obj) noexcept -> T {
    return obj.euclideanSqr();
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto atan2(const ::boyle::math::Vec2<T>& obj) noexcept -> T {
    return std::atan2(obj.y, obj.x);
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)

namespace fmt {

template <std::floating_point T, typename Char>
struct formatter<::boyle::math::Vec2<T>, Char> {
    [[using gnu: ]]
    constexpr auto parse(fmt::basic_format_parse_context<Char>& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin();
        auto end = ctx.end();
        if (it != end && *it >= '0' && *it <= '9') {
            width = 0;
            for (; it != end && *it >= '0' && *it <= '9'; ++it) {
                width = width * 10 + (*it - '0');
            }
        }
        if (it != end && *it == '.') {
            ++it;
            precision = 0;
            for (; it != end && *it >= '0' && *it <= '9'; ++it) {
                precision = precision * 10 + (*it - '0');
            }
        }
        return it;
    }

    template <typename FormatContext>
    [[using gnu: ]]
    auto format(const ::boyle::math::Vec2<T>& obj, FormatContext& ctx) const
        -> decltype(ctx.out()) {
        const std::uint8_t condition = ((width > 0) << 1) + (precision >= 0);
        std::basic_ostringstream<Char> ss;
        switch (condition) {
        case 1:
            ss << std::fixed << std::setprecision(precision) << "(x: " << obj.x << ", y: " << obj.y
               << ")";
            break;
        case 2:
            ss << std::fixed << "(x: " << std::setw(width) << obj.x << ", y: " << std::setw(width)
               << obj.y << ")";
            break;
        case 3:
            ss << std::fixed << std::setprecision(precision) << "(x: " << std::setw(width) << obj.x
               << ", y: " << std::setw(width) << obj.y << ")";
            break;
        default:
            ss << "(x: " << obj.x << ", y: " << obj.y << ")";
            break;
        }
        return fmt::format_to(ctx.out(), "{0:s}", ss.str());
    }

    int width{0};
    int precision{-1};
};

} // namespace fmt

namespace boost::serialization {

template <std::floating_point T>
[[using gnu: always_inline]]
inline constexpr auto serialize(
    auto& archive, ::boyle::math::Vec2<T>& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.x;
    archive & obj.y;
    return;
}

} // namespace boost::serialization
