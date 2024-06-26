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

    [[using gnu: always_inline]]
    Vec3() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3(T cv) noexcept
        : Triplet<T>{cv} {}
    [[using gnu: always_inline]]
    constexpr Vec3(T cx, T cy, T cz) noexcept
        : Triplet<T>{cx, cy, cz} {}
    [[using gnu: always_inline]]
    constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3(Vec3&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Vec3& other) noexcept -> Vec3& = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Vec3&& other) noexcept -> Vec3& = default;
    ~Vec3() noexcept override = default;

    template <std::floating_point U>
    [[using gnu: always_inline]]
    constexpr Vec3(const std::tuple<U>& other) noexcept
        : Triplet<T>{other} {}

    [[using gnu: always_inline]]
    constexpr Vec3(const InstanceOfTemplate<Vec3> auto& other) noexcept {
        x() = other.x();
        y() = other.y();
        z() = other.z();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto norm() const noexcept -> value_type {
        return std::hypot(x(), y(), z());
    }
    [[using gnu: pure, always_inline]]
    constexpr auto normSqr() const noexcept -> value_type {
        return x() * x() + y() * y() + z() * z();
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr auto dot(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return x() * obj.x() + y() * obj.y() + z() * obj.z();
    }
    [[using gnu: pure, always_inline]]
    constexpr auto cross(InstanceOfTemplate<Vec3> auto obj) const noexcept -> Vec3 {
        return Vec3{
            y() * obj.z() - z() * obj.y(), z() * obj.x() - x() * obj.z(),
            x() * obj.y() - y() * obj.x()
        };
    }
    [[using gnu: pure, always_inline]]
    constexpr auto euclideanTo(InstanceOfTemplate<Vec3> auto obj) const noexcept -> value_type {
        return std::hypot(x() - obj.x(), y() - obj.y(), z() - obj.z());
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-() const noexcept -> Vec3 {
        return Vec3{-x(), -y(), -z()};
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto operator==(InstanceOfTemplate<Vec3> auto other) const noexcept -> bool {
        return x() == other.x() && y() == other.y() && z() == other.z();
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator+(InstanceOfTemplate<Vec3> auto other) const noexcept -> Vec3 {
        return Vec3{x() + other.x(), y() + other.y(), z() + other.z()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator+=(InstanceOfTemplate<Vec3> auto other) noexcept -> Vec3 {
        x() += other.x();
        y() += other.y();
        z() += other.z();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator-(InstanceOfTemplate<Vec3> auto other) const noexcept -> Vec3 {
        return Vec3{x() - other.x(), y() - other.y(), z() - other.z()};
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator-=(InstanceOfTemplate<Vec3> auto other) noexcept -> Vec3& {
        x() -= other.x();
        y() -= other.y();
        z() -= other.z();
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator*(Arithmetic auto factor) const noexcept -> Vec3 {
        return Vec3{
            x() * static_cast<value_type>(factor), y() * static_cast<value_type>(factor),
            z() * static_cast<value_type>(factor)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator*=(Arithmetic auto factor) noexcept -> Vec3& {
        x() *= static_cast<value_type>(factor);
        y() *= static_cast<value_type>(factor);
        z() *= static_cast<value_type>(factor);
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto operator/(Arithmetic auto den) const noexcept -> Vec3 {
        return Vec3{
            x() / static_cast<value_type>(den), y() / static_cast<value_type>(den),
            z() / static_cast<value_type>(den)
        };
    }
    [[using gnu: always_inline, leaf]]
    constexpr auto operator/=(Arithmetic auto den) noexcept -> Vec3& {
        x() /= static_cast<value_type>(den);
        y() /= static_cast<value_type>(den);
        z() /= static_cast<value_type>(den);
        return *this;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto x() noexcept -> value_type& {
        return m_first;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto x() const noexcept -> value_type {
        return m_first;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto y() noexcept -> value_type& {
        return m_second;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto y() const noexcept -> value_type {
        return m_second;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto z() noexcept -> value_type& {
        return m_third;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr auto z() const noexcept -> value_type {
        return m_third;
    }
};

using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

[[using gnu: pure, always_inline]]
inline constexpr auto operator*(Arithmetic auto factor, InstanceOfTemplate<Vec3> auto obj) noexcept
    -> decltype(obj) {
    return obj * factor;
}

[[using gnu: pure, always_inline]]
inline constexpr auto normalize(InstanceOfTemplate<Vec3> auto obj) noexcept -> decltype(obj) {
    const typename decltype(obj)::value_type norm = obj.norm();
    return decltype(obj){obj.x() / norm, obj.y() / norm, obj.z() / norm};
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

[[using gnu: pure, always_inline]]
inline constexpr auto hypot(boyle::math::InstanceOfTemplate<boyle::math::Vec3> auto obj) noexcept ->
    typename decltype(obj)::value_type {
    return std::hypot(obj.x(), obj.y(), obj.z());
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)
