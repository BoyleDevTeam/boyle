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
#include <complex>
#include <concepts>
#include <format>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <vector>

#ifdef BOYLE_USE_AVX
#include <immintrin.h>
#endif

#include "boost/serialization/complex.hpp"
#include "fmt/format.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic T>
class Vector<T, 2> final {
  public:
    using value_type = T;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    static constexpr size_type kSize = 2;

    [[using gnu: always_inline]]
    Vector() noexcept = default;
    [[using gnu: always_inline]]
    constexpr Vector(value_type cv) noexcept
        : x{cv}, y{cv} {}
    [[using gnu: always_inline]]
    constexpr Vector(value_type cx, value_type cy) noexcept
        : x{cx}, y{cy} {}
    [[using gnu: always_inline]]
    constexpr Vector(const Vector& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(const Vector& other) noexcept -> Vector& = default;
    [[using gnu: always_inline]]
    constexpr Vector(Vector&& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr auto operator=(Vector&& other) noexcept -> Vector& = default;
    ~Vector() noexcept = default;

    template <ScalarArithmetic U>
    [[using gnu: always_inline]]
    constexpr explicit Vector(const std::pair<U, U>& other) noexcept
        : x{static_cast<value_type>(other.first)}, y{static_cast<value_type>(other.second)} {}

    template <ScalarArithmetic U>
    [[using gnu: always_inline]]
    constexpr explicit Vector(const Vector<U, 2>& other) noexcept
        : x{static_cast<value_type>(other.x)}, y{static_cast<value_type>(other.y)} {}

    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    static constexpr auto size() noexcept -> size_type {
        return kSize;
    }
    [[using gnu: const, always_inline, leaf]] [[nodiscard]]
    static constexpr auto stride() noexcept -> size_type {
        return 1;
    }

    [[using gnu: pure, always_inline]]
    constexpr auto data() noexcept -> pointer {
        return &x;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto data() const noexcept -> const_pointer {
        return &x;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_t i) noexcept -> reference {
        return (&x)[i];
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_t i) const noexcept -> const_reference {
        return (&x)[i];
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    constexpr auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= kSize) [[unlikely]] {
            throw std::out_of_range("Vector index out of range");
        }
#endif
        return (&x)[i];
    }

    [[using gnu: always_inline, hot]]
    constexpr auto updateCoeff(size_type i, const_reference value) noexcept(!BOYLE_CHECK_PARAMS)
        -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= kSize) [[unlikely]] {
            throw std::out_of_range("Vector index out of range");
        }
#endif
        return (&x)[i] = value;
    }

    [[using gnu: always_inline, leaf, hot]]
    constexpr auto operator+=(const Vector& obj) noexcept -> Vector& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int64_t> ||
                      std::is_same_v<value_type, std::uint64_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_add_epi64(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm_store_pd(
                reinterpret_cast<double*>(&x),
                _mm_add_pd(
                    _mm_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_add_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_add_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            x += obj.x;
            y += obj.y;
        }
#else
        x += obj.x;
        y += obj.y;
#endif
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator+(const Vector& obj) const& noexcept -> Vector {
        Vector result{*this};
        result += obj;
        return result;
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator+(Vector&& obj) const& noexcept -> Vector&& {
        obj += (*this);
        return std::move(obj);
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator+(const Vector& obj) && noexcept -> Vector&& {
        operator+=(obj);
        return std::move(*this);
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator+(Vector&& obj) && noexcept -> Vector&& {
        operator+=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf, hot]]
    constexpr auto operator-=(const Vector& obj) noexcept -> Vector& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int64_t> ||
                      std::is_same_v<value_type, std::uint64_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_sub_epi64(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm_store_pd(
                reinterpret_cast<double*>(&x),
                _mm_sub_pd(
                    _mm_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_sub_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_sub_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            x -= obj.x;
            y -= obj.y;
        }
#else
        x -= obj.x;
        y -= obj.y;
#endif
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator-(const Vector& obj) const& noexcept -> Vector {
        Vector result{*this};
        result -= obj;
        return result;
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator-(Vector&& obj) const& noexcept -> Vector&& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int64_t> ||
                      std::is_same_v<value_type, std::uint64_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&obj.x),
                _mm_sub_epi64(
                    _mm_load_si128(reinterpret_cast<__m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<__m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm_store_pd(
                reinterpret_cast<double*>(&obj.x),
                _mm_sub_pd(
                    _mm_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&obj.x),
                _mm_sub_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&obj.x),
                _mm256_sub_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            obj.x = x - obj.x;
            obj.y = y - obj.y;
        }
#else
        obj.x = x - obj.x;
        obj.y = y - obj.y;
#endif
        return std::move(obj);
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator-(const Vector& obj) && noexcept -> Vector&& {
        operator-=(obj);
        return std::move(*this);
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator-(Vector&& obj) && noexcept -> Vector&& {
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf, hot]]
    constexpr auto operator*=(const ScalarArithmetic auto& fac) noexcept -> Vector& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int64_t> ||
                      std::is_same_v<value_type, std::uint64_t>) {
            const auto alpha{static_cast<long long>(fac)};
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_mullo_epi64(
                    _mm_load_si128(reinterpret_cast<__m128i*>(&x)), _mm_set1_epi64x(alpha)
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            const auto alpha{static_cast<double>(fac)};
            _mm_store_pd(
                reinterpret_cast<double*>(&x),
                _mm_mul_pd(_mm_load_pd(reinterpret_cast<const double*>(&x)), _mm_set1_pd(alpha))
            );
        } else {
            const auto alpha{static_cast<value_type>(fac)};
            x *= alpha;
            y *= alpha;
        }
#else
        const auto alpha{static_cast<value_type>(fac)};
        x *= alpha;
        y *= alpha;
#endif
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator*(const ScalarArithmetic auto& fac) const& noexcept -> Vector {
        Vector result{*this};
        result *= fac;
        return result;
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator*(const ScalarArithmetic auto& fac) && noexcept -> Vector&& {
        operator*=(fac);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto operator/=(const ScalarArithmetic auto& den) noexcept -> Vector& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, double>) {
            const auto alpha{static_cast<double>(1.0 / den)};
            _mm_store_pd(
                reinterpret_cast<double*>(&x),
                _mm_mul_pd(_mm_load_pd(reinterpret_cast<const double*>(&x)), _mm_set1_pd(alpha))
            );
        } else {
            const auto alpha{static_cast<value_type>(den)};
            x /= alpha;
            y /= alpha;
        }
#else
        const auto alpha{static_cast<value_type>(den)};
        x /= alpha;
        y /= alpha;
#endif
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator/(const ScalarArithmetic auto& den) const& noexcept -> Vector {
        Vector result{*this};
        result /= den;
        return result;
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator/(const ScalarArithmetic auto& den) && noexcept -> Vector&& {
        operator/=(den);
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator-() const& noexcept -> Vector {
        return Vector{-x, -y};
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator-() && noexcept -> Vector&& {
        x = -x;
        y = -y;
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto operator==(const Vector& other) const noexcept -> bool {
        return x == other.x && y == other.y;
    }

    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto dot(const Vector& obj) const noexcept -> value_type {
        return x * obj.x + y * obj.y;
    }
    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto crossProj(const Vector& obj) const noexcept -> value_type {
        return x * obj.y - y * obj.x;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto normalized() const noexcept -> Vector {
        return *this / euclidean();
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type> {
        if constexpr (std::is_arithmetic_v<value_type>) {
            return std::hypot(x, y);
        }
        if constexpr (isComplexArithmeticV<value_type>) {
            return std::sqrt(std::norm(x) + std::norm(y));
        }
    }
    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto euclideanSqr() const noexcept -> detail::DenseNormTraitT<value_type> {
        if constexpr (std::is_arithmetic_v<value_type>) {
            return x * x + y * y;
        }
        if constexpr (isComplexArithmeticV<value_type>) {
            return std::norm(x) + std::norm(y);
        }
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclideanTo(const Vector& obj) const noexcept
        -> detail::DenseNormTraitT<value_type> {
        return operator-(obj).euclidean();
    }
    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto euclideanSqrTo(const Vector& obj) const noexcept
        -> detail::DenseNormTraitT<value_type> {
        return operator-(obj).euclideanSqr();
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto identicalTo(const Vector& obj, value_type tol = 1E-8) const noexcept -> bool {
        return euclideanSqrTo(obj) < tol * tol;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto orthogonalTo(const Vector& obj, value_type tol = 1E-8) const noexcept -> bool {
        return std::abs(dot(obj)) < tol;
    }

    [[using gnu: always_inline]]
    constexpr auto selfConjugated() noexcept -> Vector& {
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            x = std::conj(x);
            y = std::conj(y);
        }
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto conjugated() const& noexcept -> Vector {
        Matrix result{*this};
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            result.selfConjugated();
        }
        return result;
    }
    [[using gnu: always_inline]]
    constexpr auto conjugated() && noexcept -> Vector&& {
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            selfConjugated();
        }
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto angle() const noexcept -> value_type
        requires std::floating_point<value_type>
    {
        return std::atan2(y, x);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto selfRotate(const std::floating_point auto& radian) noexcept -> Vector&
        requires std::floating_point<value_type>
    {
        *this = Vector{
            x * std::cos(radian) - y * std::sin(radian), x * std::sin(radian) + y * std::cos(radian)
        };
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto rotate(const std::floating_point auto& radian) const& noexcept -> Vector
        requires std::floating_point<value_type>
    {
        return Vector{
            x * std::cos(radian) - y * std::sin(radian), x * std::sin(radian) + y * std::cos(radian)
        };
    }
    [[using gnu: always_inline, hot]]
    constexpr auto rotate(const std::floating_point auto& radian) && noexcept -> Vector&&
        requires std::floating_point<value_type>
    {
        this->selfRotate(radian);
        return std::move(*this);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto selfRotateHalfPi() noexcept -> Vector&
        requires std::floating_point<value_type>
    {
        *this = Vector{-y, x};
        return *this;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto rotateHalfPi() const& noexcept -> Vector
        requires std::floating_point<value_type>
    {
        return Vector{-y, x};
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto rotateHalfPi() && noexcept -> Vector&&
        requires std::floating_point<value_type>
    {
        this->selfRotateHalfPi();
        return std::move(*this);
    }

    alignas(2 * sizeof(value_type)) value_type x;
    value_type y;
};

template <ScalarArithmetic T>
[[using gnu: pure, always_inline, hot]]
inline constexpr auto operator*(const ScalarArithmetic auto& fac, const Vector<T, 2>& obj) noexcept
    -> Vector<T, 2> {
    return obj * fac;
}

template <ScalarArithmetic T>
[[using gnu: always_inline, hot]]
inline constexpr auto operator*(const ScalarArithmetic auto& fac, Vector<T, 2>&& obj) noexcept
    -> Vector<T, 2>&& {
    return std::move(obj) * fac;
}

template <typename Char, ScalarArithmetic T>
[[using gnu: always_inline]]
inline auto operator<<(std::basic_ostream<Char>& os, const Vector<T, 2>& obj) noexcept
    -> std::basic_ostream<Char>& {
    os << "(x: " << obj.x << ", y: " << obj.y << ")";
    return os;
}

template <ScalarArithmetic T>
using Vec2 = Vector<T, 2>;

using Vec2s = Vec2<float>;
using Vec2d = Vec2<double>;
using Vec2c = Vec2<std::complex<float>>;
using Vec2z = Vec2<std::complex<double>>;

template <typename T>
struct isVec2Arithmetic : std::false_type {};

template <ScalarArithmetic T>
struct isVec2Arithmetic<Vec2<T>> : std::true_type {};

template <typename T>
inline constexpr bool isVec2ArithmeticV = isVec2Arithmetic<T>::value;

template <typename T>
concept Vec2Arithmetic = isVec2ArithmeticV<T>;

template <Vec2Arithmetic T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline auto squeeze(
    std::ranges::forward_range auto&& xs, std::ranges::forward_range auto&& ys
) noexcept -> std::vector<T>
    requires std::same_as<std::ranges::range_value_t<decltype(xs)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(ys)>, typename T::value_type>
{
    if (xs.size() != ys.size()) [[unlikely]] {
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

template <::boyle::math::ScalarArithmetic T, typename Char>
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
inline constexpr auto hypot(const ::boyle::math::Vec2<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return std::hypot(obj.x, obj.y);
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto abs(const ::boyle::math::Vec2<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return obj.euclidean();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto norm(const ::boyle::math::Vec2<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return obj.euclideanSqr();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto conj(const ::boyle::math::Vec2<T>& vector) noexcept
    -> ::boyle::math::Vec2<T> {
    return vector.conjugated();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: always_inline]]
inline constexpr auto conj(::boyle::math::Vec2<T>&& vector) noexcept -> ::boyle::math::Vec2<T>&& {
    return std::move(vector).conjugated();
}

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto atan2(const ::boyle::math::Vec2<T>& obj) noexcept -> T {
    return std::atan2(obj.y, obj.x);
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)

namespace fmt {

template <::boyle::math::ScalarArithmetic T, typename Char>
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

template <::boyle::math::ScalarArithmetic T>
[[using gnu: always_inline]]
inline constexpr auto serialize(
    auto& archive, ::boyle::math::Vec2<T>& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.x;
    archive & obj.y;
    return;
}

} // namespace boost::serialization
