/**
 * @file vec3.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-02
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
struct Vec3 final {
    using value_type = T;

    [[using gnu: always_inline]]
    Vec3() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3(value_type cv) noexcept
        : x{cv}, y{cv}, z{cv} {}
    [[using gnu: always_inline]]
    constexpr Vec3(value_type cx, value_type cy, value_type cz) noexcept
        : x{cx}, y{cy}, z{cz} {}
    [[using gnu: always_inline]]
    constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Vec3& other) noexcept -> Vec3& = default;
    ~Vec3() noexcept = default;

    template <std::floating_point U>
    [[using gnu: always_inline]]
    constexpr explicit Vec3(const std::tuple<U>& other) noexcept
        : x{static_cast<value_type>(std::get<0>(other))},
          y{static_cast<value_type>(std::get<1>(other))},
          z{static_cast<value_type>(std::get<2>(other))} {}

    [[using gnu: always_inline]]
    constexpr Vec3(const InstanceOfTemplate<Vec3> auto& other) noexcept
        : x{static_cast<value_type>(other.x)}, y{static_cast<value_type>(other.y)},
          z{static_cast<value_type>(other.z)} {}

    [[using gnu: pure, always_inline]]
    constexpr auto euclidean() const noexcept -> value_type {
        return std::hypot(x, y, z);
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanSqr() const noexcept -> value_type {
        return x * x + y * y + z * z;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto normalized() const noexcept -> Vec3 {
        return *this / euclidean();
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return x * obj.x + y * obj.y + z * obj.z;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto cross(InstanceOfTemplate<Vec3> auto obj) const noexcept -> Vec3 {
        return Vec3{y * obj.z - z * obj.y, z * obj.x - x * obj.z, x * obj.y - y * obj.x};
    }
    [[using gnu: pure, always_inline]]
    constexpr auto crossProj(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return cross(obj).euclidean();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanTo(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return std::hypot(x - obj.x, y - obj.y, z - obj.z);
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanSqrTo(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return this->operator-(obj).euclideanSqr();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto approachTo(InstanceOfTemplate<Vec3> auto obj, value_type tol = 1E-8)
        const noexcept -> bool {
        return this->euclideanSqrTo(obj) < tol * tol;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto orthogonalTo(InstanceOfTemplate<Vec3> auto obj, value_type tol = 1E-8)
        const noexcept -> bool {
        return std::abs(this->dot(obj)) < tol;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-() const noexcept -> Vec3 {
        return Vec3{-x, -y, -z};
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto operator==(InstanceOfTemplate<Vec3> auto other) const noexcept -> bool {
        return x == other.x && y == other.y && z == other.z;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator+(InstanceOfTemplate<Vec3> auto other) const noexcept -> Vec3 {
        other.x += x;
        other.y += y;
        other.z += z;
        return other;
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator+=(InstanceOfTemplate<Vec3> auto other) noexcept -> Vec3& {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-(InstanceOfTemplate<Vec3> auto other) const noexcept -> Vec3 {
        other.x = x - other.x;
        other.y = y - other.y;
        other.z = z - other.z;
        return other;
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator-=(InstanceOfTemplate<Vec3> auto other) noexcept -> Vec3& {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator*(Arithmetic auto factor) const noexcept -> Vec3 {
        return Vec3{
            x * static_cast<value_type>(factor), y * static_cast<value_type>(factor),
            z * static_cast<value_type>(factor)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator*=(Arithmetic auto factor) noexcept -> Vec3& {
        x *= static_cast<value_type>(factor);
        y *= static_cast<value_type>(factor);
        z *= static_cast<value_type>(factor);
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator/(Arithmetic auto den) const noexcept -> Vec3 {
        return Vec3{
            x / static_cast<value_type>(den), y / static_cast<value_type>(den),
            z / static_cast<value_type>(den)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator/=(Arithmetic auto den) noexcept -> Vec3& {
        x /= static_cast<value_type>(den);
        y /= static_cast<value_type>(den);
        z /= static_cast<value_type>(den);
        return *this;
    }

    value_type x;
    value_type y;
    value_type z;
};

using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

template <typename Char>
[[using gnu: always_inline]]
inline auto operator<<(std::basic_ostream<Char>& os, InstanceOfTemplate<Vec3> auto obj) noexcept
    -> std::basic_ostream<Char>& {
    os << "(x: " << obj.x << ", y: " << obj.y << ", z: " << obj.z << ")";
    return os;
}

[[using gnu: pure, always_inline]]
inline constexpr auto operator*(Arithmetic auto factor, InstanceOfTemplate<Vec3> auto obj) noexcept
    -> decltype(obj) {
    return obj * factor;
}

template <InstanceOfTemplate<Vec3> T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline auto squeeze(
    std::ranges::forward_range auto&& xs, std::ranges::forward_range auto&& ys,
    std::ranges::forward_range auto&& zs
) noexcept -> std::vector<T>
    requires std::same_as<std::ranges::range_value_t<decltype(xs)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(ys)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(zs)>, typename T::value_type>
{
    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        return std::vector<T>{};
    }
    std::vector<T> result;
    result.reserve(xs.size());
    for (auto it_xs{xs.begin()}, it_ys{ys.begin()}, it_zs{zs.begin()};
         !(it_xs == xs.end() || it_ys == ys.end() || it_zs == zs.end());
         ++it_xs, ++it_ys, ++it_zs) {
        result.emplace_back(*it_xs, *it_ys, *it_zs);
    }
    return result;
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

template <::boyle::math::InstanceOfTemplate<::boyle::math::Vec3> T, typename Char>
struct formatter<T, Char> {
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
    auto format(T obj, FormatContext& ctx) const -> decltype(ctx.out()) {
        const std::uint8_t condition = ((width > 0) << 1) + (precision >= 0);
        std::basic_ostringstream<Char> ss;
        switch (condition) {
        case 1:
            ss << std::fixed << std::setprecision(precision) << "(x: " << obj.x << ", y: " << obj.y
               << ", z: " << obj.z << ")";
            break;
        case 2:
            ss << std::fixed << "(x: " << std::setw(width) << obj.x << ", y: " << std::setw(width)
               << obj.y << ", z: " << std::setw(width) << obj.z << ")";
            break;
        case 3:
            ss << std::fixed << std::setprecision(precision) << "(x: " << std::setw(width) << obj.x
               << ", y: " << std::setw(width) << obj.y << ", z: " << std::setw(width) << obj.z
               << ")";
            break;
        default:
            ss << "(x: " << obj.x << ", y: " << obj.y << ", z: " << obj.z << ")";
            break;
        }
        return std::format_to(ctx.out(), "{0:s}", ss.str());
    }

    int width{0};
    int precision{-1};
};

[[using gnu: pure, always_inline]]
inline constexpr auto hypot(boyle::math::InstanceOfTemplate<boyle::math::Vec3> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return std::hypot(obj.x, obj.y, obj.z);
}

[[using gnu: pure, always_inline]]
inline constexpr auto abs(boyle::math::InstanceOfTemplate<boyle::math::Vec3> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return std::hypot(obj.x, obj.y, obj.z);
}

[[using gnu: pure, always_inline]]
inline constexpr auto norm(boyle::math::InstanceOfTemplate<boyle::math::Vec3> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return obj.euclideanSqr();
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)

namespace fmt {

template <::boyle::math::InstanceOfTemplate<::boyle::math::Vec3> T, typename Char>
struct formatter<T, Char> {
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
    auto format(T obj, FormatContext& ctx) const -> decltype(ctx.out()) {
        const std::uint8_t condition = ((width > 0) << 1) + (precision >= 0);
        std::basic_ostringstream<Char> ss;
        switch (condition) {
        case 1:
            ss << std::fixed << std::setprecision(precision) << "(x: " << obj.x << ", y: " << obj.y
               << ", z: " << obj.z << ")";
            break;
        case 2:
            ss << std::fixed << "(x: " << std::setw(width) << obj.x << ", y: " << std::setw(width)
               << obj.y << ", z: " << std::setw(width) << obj.z << ")";
            break;
        case 3:
            ss << std::fixed << std::setprecision(precision) << "(x: " << std::setw(width) << obj.x
               << ", y: " << std::setw(width) << obj.y << ", z: " << std::setw(width) << obj.z
               << ")";
            break;
        default:
            ss << "(x: " << obj.x << ", y: " << obj.y << ", z: " << obj.z << ")";
            break;
        }
        return fmt::format_to(ctx.out(), "{0:s}", ss.str());
    }

    int width{0};
    int precision{-1};
};

} // namespace fmt

namespace boost::serialization {

[[using gnu: always_inline]]
inline constexpr auto serialize(
    auto& archive, boyle::math::InstanceOfTemplate<boyle::math::Vec3> auto& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.x;
    archive & obj.y;
    archive & obj.z;
    return;
}

} // namespace boost::serialization
