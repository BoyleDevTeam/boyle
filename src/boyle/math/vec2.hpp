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
    using value_type = std::remove_cv_t<T>;

    [[using gnu: always_inline]]
    Vec2() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2(T cv) noexcept
        : Duplet<T>{cv, cv} {}
    [[using gnu: always_inline]]
    constexpr Vec2(T cx, T cy) noexcept
        : Duplet<T>{cx, cy} {}
    [[using gnu: always_inline]]
    constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2(Vec2&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Vec2& other) noexcept -> Vec2& = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Vec2&& other) noexcept -> Vec2& = default;
    ~Vec2() noexcept override = default;

    template <std::floating_point U>
    [[using gnu: always_inline]]
    constexpr Vec2(const std::pair<U, U>& other) noexcept
        : Duplet<T>{other} {}

    [[using gnu: always_inline]]
    constexpr Vec2(const InstanceOfTemplate<Vec2> auto& other) noexcept {
        x() = other.x();
        y() = other.y();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto norm() const noexcept -> value_type {
        return std::hypot(x(), y());
    }
    [[using gnu: pure, always_inline]]
    constexpr auto normSqr() const noexcept -> value_type {
        return x() * x() + y() * y();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto angle() const noexcept -> value_type {
        return std::atan2(y(), x());
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(InstanceOfTemplate<Vec2> auto obj) const noexcept -> value_type {
        return x() * obj.x() + y() * obj.y();
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto cross(InstanceOfTemplate<Vec2> auto obj) const noexcept -> value_type {
        return x() * obj.y() - y() * obj.x();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanTo(InstanceOfTemplate<Vec2> auto obj) const noexcept -> value_type {
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
    constexpr auto operator==(InstanceOfTemplate<Vec2> auto other) const noexcept -> bool {
        return x() == other.x() && y() == other.y();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-() const noexcept -> Vec2 {
        return Vec2{-x(), -y()};
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator+(InstanceOfTemplate<Vec2> auto other) const noexcept -> Vec2 {
        return Vec2{x() + other.x(), y() + other.y()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator+=(InstanceOfTemplate<Vec2> auto other) noexcept -> Vec2& {
        x() += other.x();
        y() += other.y();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-(InstanceOfTemplate<Vec2> auto other) const noexcept -> Vec2 {
        return Vec2{x() - other.x(), y() - other.y()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator-=(InstanceOfTemplate<Vec2> auto other) noexcept -> Vec2& {
        x() -= other.x();
        y() -= other.y();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator*(Arithmetic auto factor) const noexcept -> Vec2 {
        return Vec2{x() * static_cast<value_type>(factor), y() * static_cast<value_type>(factor)};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator*=(Arithmetic auto factor) noexcept -> Vec2& {
        x() *= static_cast<value_type>(factor);
        y() *= static_cast<value_type>(factor);
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator/(Arithmetic auto den) const noexcept -> Vec2 {
        return Vec2{x() / static_cast<value_type>(den), y() / static_cast<value_type>(den)};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator/=(Arithmetic auto den) noexcept -> Vec2& {
        x() /= static_cast<value_type>(den);
        y() /= static_cast<value_type>(den);
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

[[using gnu: pure, always_inline]]
inline constexpr auto operator*(Arithmetic auto factor, InstanceOfTemplate<Vec2> auto obj) noexcept
    -> decltype(obj) {
    return obj * factor;
}

[[using gnu: pure, always_inline]]
inline constexpr auto normalize(InstanceOfTemplate<Vec2> auto obj) noexcept -> decltype(obj) {
    const typename decltype(obj)::value_type norm = obj.norm();
    return decltype(obj){obj.x() / norm, obj.y() / norm};
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

[[using gnu: pure, always_inline]]
inline constexpr auto hypot(boyle::math::InstanceOfTemplate<boyle::math::Vec2> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return std::hypot(obj.x(), obj.y());
}

[[using gnu: pure, always_inline]]
inline constexpr auto atan2(boyle::math::InstanceOfTemplate<boyle::math::Vec2> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return std::atan2(obj.y(), obj.x());
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
