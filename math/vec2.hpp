/**
 * @file vec2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-17
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

enum class Vec2Mode : unsigned int {
    XY,
    SL,
    ST,
    LT
};

template <Vec2Mode VM, typename T>
class Vec2;

template <typename T>
class Vec2<Vec2Mode::XY, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : x(cv), y(cv) {}
    [[using gnu: always_inline]] constexpr Vec2(T cx, T cy) noexcept : x(cx), y(cy) {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2&
    operator=(const Vec2& other) noexcept = default;
    ~Vec2() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec2<Vec2Mode::XY, U>() noexcept {
        return Vec2<Vec2Mode::XY, U>(static_cast<U>(x), static_cast<U>(y));
    }

    [[using gnu: pure, always_inline]]
    constexpr T length() const noexcept {
        return std::hypot(x, y);
    }
    [[using gnu: pure, always_inline]]
    constexpr T lengthSqr() const noexcept {
        return x * x + y * y;
    }
    [[using gnu: pure, always_inline]]
    constexpr T angle() const noexcept {
        return std::atan2(y, x);
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T dot(Vec2 obj) const noexcept {
        return x * obj.x + y * obj.y;
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T cross(Vec2 obj) const noexcept {
        return x * obj.y - y * obj.x;
    }
    [[using gnu: pure, always_inline]]
    constexpr T distanceTo(Vec2 obj) const noexcept {
        return std::hypot(x - obj.x, y - obj.y);
    }

    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2 rotate(U radian) const noexcept {
        static_assert(
            std::is_floating_point_v<U>, "The loaded type must be a floating-point type."
        );
        return Vec2{
            x * std::cos(radian) - y * std::sin(radian),
            x * std::sin(radian) + y * std::cos(radian)};
    }
    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2& selfRotate(U radian) noexcept {
        static_assert(
            std::is_floating_point_v<U>, "The loaded type must be a floating-point type."
        );
        *this = Vec2{
            x * std::cos(radian) - y * std::sin(radian),
            x * std::sin(radian) + y * std::cos(radian)};
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2 rotateHalfPi() const noexcept {
        return Vec2{-y, x};
    }
    [[using gnu: always_inline]]
    constexpr Vec2& selfRotateHalfPi() noexcept {
        *this = Vec2{-y, x};
        return *this;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr bool
    operator==(Vec2 other) {
        return x == other.x && y == other.y;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator-() const noexcept {
        return Vec2{-x, -y};
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator+(Vec2 other) const noexcept {
        return Vec2{x + other.x, y + other.y};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator+=(Vec2 other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator-(Vec2 other) const noexcept {
        return Vec2{x - other.x, y - other.y};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator-=(Vec2 other) noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator*(U factor) const noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        return Vec2{x * static_cast<T>(factor), y * static_cast<T>(factor)};
    }
    template <typename U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator*=(U factor) noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        x *= static_cast<T>(factor);
        y *= static_cast<T>(factor);
        return *this;
    }

    template <typename U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator/(U den) const noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        return Vec2{x / static_cast<T>(den), y / static_cast<T>(den)};
    }
    template <typename U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator/=(U den) noexcept {
        static_assert(std::is_arithmetic_v<U>, "The loaded type must has arithmetic operators.");
        x /= static_cast<T>(den);
        y /= static_cast<T>(den);
        return *this;
    }

    T x;
    T y;
};

template <typename T>
class Vec2<Vec2Mode::SL, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : s(cv), l(cv) {}
    [[using gnu: always_inline]] constexpr Vec2(T cs, T cl) noexcept : s(cs), l(cl) {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2&
    operator=(const Vec2& other) noexcept = default;
    ~Vec2() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec2<Vec2Mode::SL, U>() noexcept {
        return Vec2<Vec2Mode::SL, U>{static_cast<U>(s), static_cast<U>(l)};
    }

    T s;
    T l;
};

template <typename T>
class Vec2<Vec2Mode::ST, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : s(cv), t(cv) {}
    [[using gnu: always_inline]] constexpr Vec2(T cs, T ct) noexcept : s(cs), t(ct) {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2&
    operator=(const Vec2& other) noexcept = default;
    ~Vec2() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec2<Vec2Mode::SL, U>() noexcept {
        return Vec2<Vec2Mode::ST, U>{static_cast<U>(s), static_cast<U>(t)};
    }

    T s;
    T t;
};

template <typename T>
class Vec2<Vec2Mode::LT, T> final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept = default;
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : l(cv), t(cv) {}
    [[using gnu: always_inline]] constexpr Vec2(T cl, T ct) noexcept : l(cl), t(ct) {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vec2&
    operator=(const Vec2& other) noexcept = default;
    ~Vec2() noexcept = default;

    template <typename U>
    [[using gnu: pure, always_inline]] constexpr operator Vec2<Vec2Mode::LT, U>() noexcept {
        return Vec2<Vec2Mode::ST, U>{static_cast<U>(l), static_cast<U>(t)};
    }

    T l;
    T t;
};

using Vec2f = Vec2<Vec2Mode::XY, float>;
using Vec2d = Vec2<Vec2Mode::XY, double>;

using SlVec2f = Vec2<Vec2Mode::SL, float>;
using SlVec2d = Vec2<Vec2Mode::SL, double>;

using StVec2f = Vec2<Vec2Mode::ST, float>;
using StVec2d = Vec2<Vec2Mode::ST, double>;

using LtVec2f = Vec2<Vec2Mode::LT, float>;
using LtVec2d = Vec2<Vec2Mode::LT, double>;

template <typename T, typename U>
[[using gnu: const, always_inline]]
constexpr inline Vec2<Vec2Mode::XY, T>
operator*(U factor, Vec2<Vec2Mode::XY, T> obj) noexcept {
    return obj * factor;
}

template <typename T>
[[using gnu: const, always_inline]]
constexpr inline Vec2<Vec2Mode::XY, T> normalize(Vec2<Vec2Mode::XY, T> obj) noexcept {
    const T length = obj.length();
    return Vec2<Vec2Mode::XY, T>{obj.x / length, obj.y / length};
}

template <typename T>
struct isVecArithmetic<tiny_pnc::math::Vec2<Vec2Mode::XY, T>> {
    static constexpr bool value = true;
};

} // namespace math
} // namespace tiny_pnc

namespace std {

template <typename T>
[[using gnu: const, always_inline]]
constexpr inline T hypot(tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> obj) noexcept {
    return std::hypot(obj.x, obj.y);
}

template <typename T>
[[using gnu: const, always_inline]]
constexpr inline T atan2(tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> obj) noexcept {
    return std::atan2(obj.y, obj.x);
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T> obj) noexcept {
    os << "(" << obj.x << ", " << obj.y << ")";
    return os;
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T> obj) noexcept {
    os << "(" << obj.s << ", " << obj.l << ")";
    return os;
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::ST, T> obj) noexcept {
    os << "(" << obj.s << ", " << obj.t << ")";
    return os;
}

template <typename T>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::LT, T> obj) noexcept {
    os << "(" << obj.l << ", " << obj.t << ")";
    return os;
}

} // namespace std

namespace boost {
namespace serialization {

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::XY, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.x;
    ar& obj.y;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::SL, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.s;
    ar& obj.l;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::ST, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.s;
    ar& obj.t;
    return;
}

template <typename Archive, typename T>
[[using gnu: always_inline]]
constexpr inline void serialize(
    Archive& ar, tiny_pnc::math::Vec2<tiny_pnc::math::Vec2Mode::LT, T>& obj,
    const unsigned int version
) noexcept {
    ar& obj.l;
    ar& obj.t;
    return;
}

} // namespace serialization
} // namespace boost
