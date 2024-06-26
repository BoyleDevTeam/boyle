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

#include "boyle/math/concepts.hpp"
#include "boyle/math/triplet.hpp"

namespace boyle::math {

template <std::floating_point T>
class Vec3 final : public Triplet<T> {
    template <unsigned int Num, typename TripletType>
    friend struct TripletElement;
    using Triplet<T>::m_first;
    using Triplet<T>::m_second;
    using Triplet<T>::m_third;

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec3() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec3(T cv) noexcept : Triplet<T>{cv} {}
    [[using gnu: always_inline]] constexpr Vec3(T cx, T cy, T cz) noexcept
        : Triplet<T>{cx, cy, cz} {}
    [[using gnu: always_inline]] constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]] constexpr Vec3(Vec3&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto
    operator=(const Vec3& other) noexcept -> Vec3& = default;
    [[using gnu: always_inline]]
    constexpr auto
    operator=(Vec3&& other) noexcept -> Vec3& = default;
    ~Vec3() noexcept override = default;

    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr
    operator Vec3<U>() noexcept {
        return Vec3<U>{static_cast<U>(x()), static_cast<U>(y()), static_cast<U>(z())};
    }

    [[using gnu: pure, always_inline]]
    constexpr auto length() const noexcept -> T {
        return std::hypot(x(), y(), z());
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(Vec3 obj) const noexcept -> T {
        return x() * obj.x() + y() * obj.y() + z() * obj.z();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto cross(Vec3 obj) const noexcept -> Vec3 {
        return Vec3{
            y() * obj.z() - z() * obj.y(), z() * obj.x() - x() * obj.z(),
            x() * obj.y() - y() * obj.x()
        };
    }
    [[using gnu: pure, always_inline]]
    constexpr auto distanceTo(Vec3 obj) const noexcept -> T {
        return std::hypot(x() - obj.x(), y() - obj.y(), z() - obj.z());
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator-() const noexcept -> Vec3 {
        return Vec3{-x(), -y(), -z()};
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto
    operator==(Vec3 other) const noexcept -> bool {
        return x() == other.x() && y() == other.y() && z() == other.z();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator+(Vec3 other) const noexcept -> Vec3 {
        return Vec3{x() + other.x(), y() + other.y(), z() + other.z()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator+=(Vec3 other) noexcept -> Vec3 {
        x() += other.x();
        y() += other.y();
        z() += other.z();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator-(Vec3 other) const noexcept -> Vec3 {
        return Vec3{x() - other.x(), y() - other.y(), z() - other.z()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator-=(Vec3 other) noexcept -> Vec3& {
        x() -= other.x();
        y() -= other.y();
        z() -= other.z();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator*(Arithmetic auto factor) const noexcept -> Vec3 {
        return Vec3{
            x() * static_cast<T>(factor), y() * static_cast<T>(factor), z() * static_cast<T>(factor)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator*=(Arithmetic auto factor) noexcept -> Vec3& {
        x() *= static_cast<T>(factor);
        y() *= static_cast<T>(factor);
        z() *= static_cast<T>(factor);
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto
    operator/(Arithmetic auto den) const noexcept -> Vec3 {
        return Vec3{
            x() / static_cast<T>(den), y() / static_cast<T>(den), z() / static_cast<T>(den)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto
    operator/=(Arithmetic auto den) noexcept -> Vec3& {
        x() /= static_cast<T>(den);
        y() /= static_cast<T>(den);
        z() /= static_cast<T>(den);
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

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto z() noexcept -> T& {
        return m_third;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto z() const noexcept -> T {
        return m_third;
    }
};

using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto
operator*(Arithmetic auto factor, Vec3<T> obj) noexcept -> Vec3<T> {
    return obj * factor;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto normalize(Vec3<T> obj) noexcept -> Vec3<T> {
    T length = obj.length();
    return Vec3<T>{obj.x() / length, obj.y() / length, obj.z() / length};
}

} // namespace boyle::math

namespace std {

template <floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto hypot(boyle::math::Vec3<T> obj) noexcept -> T {
    return std::hypot(obj.x(), obj.y(), obj.z());
}

} // namespace std
