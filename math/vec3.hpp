/**
 * @file vec3.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-02
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <cmath>
#include <ostream>
#include <type_traits>

#include "math/type_traits.hpp"

namespace tiny_pnc {
namespace math {

enum class Vec3Mode : unsigned int {
    XYZ,
    SLZ,
    SLT
};

template <Vec3Mode VM, typename T>
class Vec3;

template <typename T>
class Vec3<Vec3Mode::XYZ, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec3() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec3(T cv) noexcept : x(cv), y(cv), z(cv) {}
    [[using gnu: always_inline]] constexpr Vec3(T cx, T cy, T cz) noexcept : x(cx), y(cy), z(cz) {}
    [[using gnu: always_inline]] constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3&
    operator=(const Vec3& other) noexcept = default;
    ~Vec3() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec3<Vec3Mode::XYZ, U>() noexcept {
        return Vec3<Vec3Mode::XYZ, U>{static_cast<U>(x), static_cast<U>(y), static_cast<U>(z)};
    }

    [[using gnu: pure, always_inline]]
    constexpr T length() const noexcept {
        return std::hypot(x, y, z);
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T dot(Vec3 obj) const noexcept {
        return x * obj.x + y * obj.y + z * obj.z;
    }
    [[using gnu: pure, always_inline]]
    constexpr Vec3 cross(Vec3 obj) const noexcept {
        return Vec3{y * obj.z - z * obj.y, z * obj.x - x * obj.z, x * obj.y - y * obj.x};
    }
    [[using gnu: pure, always_inline]]
    constexpr T distanceTo(Vec3 obj) const noexcept {
        return std::hypot(x - obj.x, y - obj.y, z - obj.z);
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator-() const noexcept {
        return Vec3{-x, -y, -z};
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr bool
    operator==(Vec3 other) {
        return x == other.x && y == other.y && z == other.z;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator+(Vec3 other) const noexcept {
        return Vec3{x + other.x, y + other.y, z + other.z};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator+=(Vec3 other) noexcept {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator-(Vec3 other) const noexcept {
        return Vec3{x - other.x, y - other.y, z - other.z};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator-=(Vec3 other) noexcept {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator*(U factor) const noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        return Vec3{
            x * static_cast<T>(factor), y * static_cast<T>(factor), z * static_cast<T>(factor)};
    }
    template <typename U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator*=(U factor) noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        x *= static_cast<T>(factor);
        y *= static_cast<T>(factor);
        z *= static_cast<T>(factor);
        return *this;
    }

    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator/(U den) const noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        return Vec3{x / static_cast<T>(den), y / static_cast<T>(den), z / static_cast<T>(den)};
    }
    template <typename U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator/=(U den) noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        x /= static_cast<T>(den);
        y /= static_cast<T>(den);
        z /= static_cast<T>(den);
        return *this;
    }

    T x;
    T y;
    T z;
};

template <typename T>
class Vec3<Vec3Mode::SLZ, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec3() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec3(T cv) noexcept : s(cv), l(cv), z(cv) {}
    [[using gnu: always_inline]] constexpr Vec3(T cs, T cl, T cz) noexcept : s(cs), l(cl), z(cz) {}
    [[using gnu: always_inline]] constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3&
    operator=(const Vec3& other) noexcept = default;
    ~Vec3() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec3<Vec3Mode::SLZ, U>() noexcept {
        return Vec3<Vec3Mode::SLZ, U>{static_cast<U>(s), static_cast<U>(l), static_cast<U>(z)};
    }

    T s;
    T l;
    T z;
};

template <typename T>
class Vec3<Vec3Mode::SLT, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec3() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec3(T cv) noexcept : s(cv), l(cv), t(cv) {}
    [[using gnu: always_inline]] constexpr Vec3(T cs, T cl, T ct) noexcept : s(cs), l(cl), t(ct) {}
    [[using gnu: always_inline]] constexpr Vec3(const Vec3& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec3&
    operator=(const Vec3& other) noexcept = default;
    ~Vec3() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec3<Vec3Mode::SLT, U>() noexcept {
        return Vec3<Vec3Mode::SLT, U>{static_cast<U>(s), static_cast<U>(l), static_cast<U>(t)};
    }

    T s;
    T l;
    T t;
};

using Vec3f = Vec3<Vec3Mode::XYZ, float>;
using Vec3d = Vec3<Vec3Mode::XYZ, double>;

using SlzVec3f = Vec3<Vec3Mode::SLZ, float>;
using SlzVec3d = Vec3<Vec3Mode::SLZ, double>;

using SltVec3f = Vec3<Vec3Mode::SLT, float>;
using SltVec3d = Vec3<Vec3Mode::SLT, double>;

template <typename T, typename U>
[[using gnu: const, always_inline]]
constexpr inline Vec3<Vec3Mode::XYZ, T>
operator*(U factor, Vec3<Vec3Mode::XYZ, T> obj) noexcept {
    return obj * factor;
}

template <typename T>
[[using gnu: const, always_inline]]
constexpr inline Vec3<Vec3Mode::XYZ, T> normalize(Vec3<Vec3Mode::XYZ, T> obj) noexcept {
    T length = obj.length();
    return Vec3<Vec3Mode::XYZ, T>{obj.x / length, obj.y / length, obj.z / length};
}

template <typename T>
struct isVecArithmetic<tiny_pnc::math::Vec3<Vec3Mode::XYZ, T>> {
    static constexpr bool value = true;
};

} // namespace math
} // namespace tiny_pnc

namespace std {

template <typename T>
[[using gnu: pure, always_inline]]
constexpr inline T hypot(tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::XYZ, T> obj) noexcept {
    return std::hypot(obj.x, obj.y, obj.z);
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::XYZ, T> obj) noexcept {
    os << "(" << obj.x << ", " << obj.y << ", " << obj.z << ")";
    return os;
}

template <typename T>
[[using gnu: pure, always_inline]]
constexpr inline T hypot(tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLZ, T> obj) noexcept {
    return std::hypot(obj.s, obj.l, obj.z);
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLZ, T> obj) noexcept {
    os << "(" << obj.s << ", " << obj.l << ", " << obj.z << ")";
    return os;
}

template <typename T>
[[using gnu: pure, always_inline]]
constexpr inline T hypot(tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLT, T> obj) noexcept {
    return std::hypot(obj.s, obj.l, obj.t);
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLT, T> obj) noexcept {
    os << "(" << obj.s << ", " << obj.l << ", " << obj.t << ")";
    return os;
}

} // namespace std

namespace boost {
namespace serialization {

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::XYZ, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.x;
    ar& obj.y;
    ar& obj.z;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLZ, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.s;
    ar& obj.l;
    ar& obj.z;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec3<tiny_pnc::math::Vec3Mode::SLT, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.s;
    ar& obj.l;
    ar& obj.t;
    return;
}

} // namespace serialization
} // namespace boost
