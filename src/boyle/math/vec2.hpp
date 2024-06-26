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

#include "boyle/math/concepts.hpp"
#include "boyle/math/duplet.hpp"

namespace boyle::math {

template <std::floating_point T>
class Vec2 final : public Duplet<T> {
    template <unsigned int Num, typename DupletType>
    friend struct DupletElement;
    using Duplet<T>::m_first;
    using Duplet<T>::m_second;

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : Duplet<T>{cv, cv} {}
    [[using gnu: always_inline]] constexpr Vec2(T cx, T cy) noexcept : Duplet<T>{cx, cy} {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(Vec2&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto
    operator=(const Vec2& other) noexcept -> Vec2& = default;
    [[using gnu: always_inline]]
    constexpr auto
    operator=(Vec2&& other) noexcept -> Vec2& = default;
    ~Vec2() noexcept override = default;

    template <std::floating_point U>
    [[using gnu: always_inline]] constexpr Vec2(const std::pair<U, U>& other) noexcept
        : Duplet<T>{other} {}

    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr
    operator Vec2<U>() const noexcept {
        return Vec2<U>{static_cast<U>(x()), static_cast<U>(y())};
    }

    [[using gnu: pure, always_inline]]
    constexpr auto length() const noexcept -> T {
        return std::hypot(x(), y());
    }
    [[using gnu: pure, always_inline]]
    constexpr auto lengthSqr() const noexcept -> T {
        return x() * x() + y() * y();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto angle() const noexcept -> T {
        return std::atan2(y(), x());
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(Vec2 obj) const noexcept -> T {
        return x() * obj.x() + y() * obj.y();
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto cross(Vec2 obj) const noexcept -> T {
        return x() * obj.y() - y() * obj.x();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto distanceTo(Vec2 obj) const noexcept -> T {
        return std::hypot(x() - obj.x(), y() - obj.y());
    }

    [[using gnu: pure, always_inline]]
    constexpr auto rotate(std::floating_point auto radian) const noexcept -> Vec2 {
        return Vec2{
            x() * std::cos(radian) - y() * std::sin(radian),
            x() * std::sin(radian) + y() * std::cos(radian)
        };
    }
    [[using gnu: pure, always_inline]]
    constexpr auto selfRotate(std::floating_point auto radian) noexcept -> Vec2& {
        *this = Vec2{
            x() * std::cos(radian) - y() * std::sin(radian),
            x() * std::sin(radian) + y() * std::cos(radian)
        };
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto rotateHalfPi() const noexcept -> Vec2 {
        return Vec2{-y(), x()};
    }
    [[using gnu: always_inline]]
    constexpr auto selfRotateHalfPi() noexcept -> Vec2& {
        *this = Vec2{-y(), x()};
        return *this;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto
    operator==(Vec2 other) const noexcept -> bool {
        return x() == other.x() && y() == other.y();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator-() const noexcept -> Vec2 {
        return Vec2{-x(), -y()};
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator+(Vec2 other) const noexcept -> Vec2 {
        return Vec2{x() + other.x(), y() + other.y()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator+=(Vec2 other) noexcept -> Vec2& {
        x() += other.x();
        y() += other.y();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator-(Vec2 other) const noexcept -> Vec2 {
        return Vec2{x() - other.x(), y() - other.y()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator-=(Vec2 other) noexcept -> Vec2& {
        x() -= other.x();
        y() -= other.y();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator*(Arithmetic auto factor) const noexcept -> Vec2 {
        return Vec2{x() * static_cast<T>(factor), y() * static_cast<T>(factor)};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator*=(Arithmetic auto factor) noexcept -> Vec2& {
        x() *= static_cast<T>(factor);
        y() *= static_cast<T>(factor);
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator/(Arithmetic auto den) const noexcept -> Vec2 {
        return Vec2{x() / static_cast<T>(den), y() / static_cast<T>(den)};
    }
    template <Arithmetic U>
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator/=(Arithmetic auto den) noexcept -> Vec2& {
        x() /= static_cast<T>(den);
        y() /= static_cast<T>(den);
        return *this;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto x() noexcept -> T& {
        return m_first;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto x() const noexcept -> T {
        return m_first;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto y() noexcept -> T& {
        return m_second;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto y() const noexcept -> T {
        return m_second;
    }
};

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto
operator*(Arithmetic auto factor, Vec2<T> obj) noexcept -> Vec2<T> {
    return obj * factor;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto normalize(Vec2<T> obj) noexcept -> Vec2<T> {
    const T length = obj.length();
    return Vec2<T>{obj.x() / length, obj.y() / length};
}

} // namespace boyle::math

namespace std {

template <floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto hypot(boyle::math::Vec2<T> obj) noexcept -> T {
    return std::hypot(obj.x(), obj.y());
}

template <floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto atan2(boyle::math::Vec2<T> obj) noexcept -> T {
    return std::atan2(obj.y(), obj.x());
}

} // namespace std
