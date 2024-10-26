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
class alignas(4 * sizeof(T)) Vector<T, 3> final {
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
        : x{cv}, y{cv}, z{cv} {}
    [[using gnu: always_inline]]
    constexpr Vector(value_type cx, value_type cy, value_type cz) noexcept
        : x{cx}, y{cy}, z{cz} {}
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
    constexpr explicit Vector(const std::tuple<U, U, U>& other) noexcept
        : x{static_cast<value_type>(std::get<0>(other))},
          y{static_cast<value_type>(std::get<1>(other))},
          z{static_cast<value_type>(std::get<2>(other))} {}

    template <ScalarArithmetic U>
    [[using gnu: always_inline]]
    constexpr explicit Vector(const Vector<U, 3>& other) noexcept
        : x{static_cast<value_type>(other.x)}, y{static_cast<value_type>(other.y)},
          z{static_cast<value_type>(other.z)} {}

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
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_add_epi32(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&x),
                _mm256_add_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)),
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_add_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_add_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm256_store_ps(
                reinterpret_cast<float*>(&x),
                _mm256_add_ps(
                    _mm256_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm256_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm512_store_pd(
                reinterpret_cast<double*>(&x),
                _mm512_add_pd(
                    _mm512_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm512_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            x += obj.x;
            y += obj.y;
            z += obj.z;
        }
#else
        x += obj.x;
        y += obj.y;
        z += obj.z;
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
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_sub_epi32(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&x),
                _mm256_sub_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)),
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_sub_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_sub_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm256_store_ps(
                reinterpret_cast<float*>(&x),
                _mm256_sub_ps(
                    _mm256_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm256_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm512_store_pd(
                reinterpret_cast<double*>(&x),
                _mm512_sub_pd(
                    _mm512_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm512_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            x -= obj.x;
            y -= obj.y;
            z -= obj.z;
        }
#else
        x -= obj.x;
        y -= obj.y;
        z -= obj.z;
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
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&obj.x),
                _mm_sub_epi32(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&obj.x),
                _mm256_sub_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)),
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&obj.x),
                _mm_sub_ps(
                    _mm_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&obj.x),
                _mm256_sub_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm256_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm256_store_ps(
                reinterpret_cast<float*>(&obj.x),
                _mm256_sub_ps(
                    _mm256_load_ps(reinterpret_cast<const float*>(&x)),
                    _mm256_load_ps(reinterpret_cast<const float*>(&obj.x))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm512_store_pd(
                reinterpret_cast<double*>(&obj.x),
                _mm512_sub_pd(
                    _mm512_load_pd(reinterpret_cast<const double*>(&x)),
                    _mm512_load_pd(reinterpret_cast<const double*>(&obj.x))
                )
            );
        } else {
            obj.x = x - obj.x;
            obj.y = y - obj.y;
            obj.z = z - obj.z;
        }
#else
        obj.x = x - obj.x;
        obj.y = y - obj.y;
        obj.z = z - obj.z;
#endif
        return std::move(obj);
    }
    [[using gnu: always_inline, hot]]
    constexpr auto operator-(const Vector& obj) && noexcept -> Vector&& {
        operator-=(obj);
        return std::move(*this);
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator-(Vector&& obj) && noexcept -> Vector&& {
        operator-=(obj);
        return std::move(*this);
    }

    [[using gnu: always_inline, leaf, hot]]
    constexpr auto operator*=(const ScalarArithmetic auto& fac) noexcept -> Vector& {
        const auto alpha{static_cast<value_type>(fac)};
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_mullo_epi32(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)), _mm_set1_epi32(alpha)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&x),
                _mm256_mullo_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)),
                    _mm256_set1_epi64x(alpha)
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_mul_ps(_mm_load_ps(reinterpret_cast<const float*>(&x)), _mm_set1_ps(alpha))
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_mul_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)), _mm256_set1_pd(alpha)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            const __m256 lhs{_mm256_load_ps(reinterpret_cast<const float*>(&x))};
            _mm256_store_ps(
                reinterpret_cast<float*>(&x),
                _mm256_addsub_ps(
                    _mm256_mul_ps(lhs, _mm256_set1_ps(alpha.real())),
                    _mm256_mul_ps(_mm256_permute_ps(lhs, 0xB1), _mm256_set1_ps(alpha.imag()))
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            const __m512d lhs{_mm512_load_pd(reinterpret_cast<const double*>(&x))};
            _mm512_store_pd(
                reinterpret_cast<double*>(&x),
                _mm512_add_pd(
                    _mm512_mul_pd(lhs, _mm512_set1_pd(alpha.real())),
                    _mm512_xor_pd(
                        _mm512_mul_pd(_mm512_permute_pd(lhs, 0x55), _mm512_set1_pd(alpha.imag())),
                        _mm512_set_pd(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0)
                    )
                )
            );
        } else {
            x *= alpha;
            y *= alpha;
            z *= alpha;
        }
#else
        x *= alpha;
        y *= alpha;
        z *= alpha;
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
        if constexpr (std::is_integral_v<value_type>) {
            const auto alpha{static_cast<value_type>(den)};
            x /= alpha;
            y /= alpha;
            z /= alpha;
        } else {
            const auto alpha{static_cast<value_type>(1.0) / static_cast<value_type>(den)};
            operator*=(alpha);
        }
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
        Vector result{*this};
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&result.x),
                _mm_xor_epi32(
                    _mm_load_si128(reinterpret_cast<const __m128i*>(&x)), _mm_set1_epi32(-0)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&result.x),
                _mm256_xor_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)), _mm256_set1_epi64x(-0L)
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&result.x),
                _mm_xor_ps(_mm_load_ps(reinterpret_cast<const float*>(&x)), _mm_set1_ps(-0.0F))
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&result.x),
                _mm256_xor_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)), _mm256_set1_pd(-0.0)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm256_store_ps(
                reinterpret_cast<float*>(&result.x),
                _mm256_xor_ps(
                    _mm256_load_ps(reinterpret_cast<const float*>(&x)), _mm256_set1_ps(-0.0F)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm512_store_pd(
                reinterpret_cast<double*>(&result.x),
                _mm512_xor_pd(
                    _mm512_load_pd(reinterpret_cast<const double*>(&x)), _mm512_set1_pd(-0.0)
                )
            );
        } else {
            result.x = -x;
            result.y = -y;
            result.z = -z;
        }
#else
        result.x = -x;
        result.y = -y;
        result.z = -z;
#endif
        return result;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator-() && noexcept -> Vector&& {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&x),
                _mm_xor_epi32(_mm_load_si128(reinterpret_cast<__m128i*>(&x)), _mm_set1_epi32(-0))
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&x),
                _mm256_xor_epi64(
                    _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)), _mm256_set1_epi64x(-0L)
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            _mm_store_ps(
                reinterpret_cast<float*>(&x),
                _mm_xor_ps(_mm_load_ps(reinterpret_cast<const float*>(&x)), _mm_set1_ps(-0.0F))
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            _mm256_store_pd(
                reinterpret_cast<double*>(&x),
                _mm256_xor_pd(
                    _mm256_load_pd(reinterpret_cast<const double*>(&x)), _mm256_set1_pd(-0.0)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            _mm256_store_ps(
                reinterpret_cast<float*>(&x),
                _mm256_xor_ps(
                    _mm256_load_ps(reinterpret_cast<const float*>(&x)), _mm256_set1_ps(-0.0F)
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            _mm512_store_pd(
                reinterpret_cast<double*>(&x),
                _mm512_xor_pd(
                    _mm512_load_pd(reinterpret_cast<const double*>(&x)), _mm512_set1_pd(-0.0)
                )
            );
        } else {
            x = -x;
            y = -y;
            z = -z;
        }
#else
        x = -x;
        y = -y;
        z = -z;
#endif
        return std::move(*this);
    }

    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto operator==(const Vector& other) const noexcept -> bool {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            return _mm_cmp_epi32_mask(
                       _mm_load_si128(reinterpret_cast<const __m128i*>(&x)),
                       _mm_load_si128(reinterpret_cast<const __m128i*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x7;
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            return _mm256_cmp_epi64_mask(
                       _mm256_load_si256(reinterpret_cast<const __m256i*>(&x)),
                       _mm256_load_si256(reinterpret_cast<const __m256i*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x7;
        } else if constexpr (std::is_same_v<value_type, float>) {
            return _mm_cmp_ps_mask(
                       _mm_load_ps(reinterpret_cast<const float*>(&x)),
                       _mm_load_ps(reinterpret_cast<const float*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x7;
        } else if constexpr (std::is_same_v<value_type, double>) {
            return _mm256_cmp_pd_mask(
                       _mm256_load_pd(reinterpret_cast<const double*>(&x)),
                       _mm256_load_pd(reinterpret_cast<const double*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x7;
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            return _mm256_cmp_ps_mask(
                       _mm256_load_ps(reinterpret_cast<const float*>(&x)),
                       _mm256_load_ps(reinterpret_cast<const float*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x3F;
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            return _mm512_cmp_pd_mask(
                       _mm512_load_pd(reinterpret_cast<const double*>(&x)),
                       _mm512_load_pd(reinterpret_cast<const double*>(&other.x)), _CMP_EQ_OQ
                   ) == 0x3F;
        } else {
            return x == other.x && y == other.y && z == other.z;
        }
#else
        return x == other.x && y == other.y && z == other.z;
#endif
    }

    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto dot(const Vector& obj) const noexcept -> value_type {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            __m128i sum{_mm_mullo_epi32(
                _mm_maskz_load_epi32(0b0111, reinterpret_cast<const std::int32_t*>(&x)),
                _mm_maskz_load_epi32(0b0111, reinterpret_cast<const std::int32_t*>(&obj.x))
            )};
            sum = _mm_hadd_epi32(sum, sum);
            return _mm_cvtsi128_si32(_mm_hadd_epi32(sum, sum));
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            const __m256i sum{_mm256_mullo_epi64(
                _mm256_maskz_load_epi64(0b0111, reinterpret_cast<const std::int64_t*>(&x)),
                _mm256_maskz_load_epi64(0b0111, reinterpret_cast<const std::int64_t*>(&obj.x))
            )};
            const __m128i sum2{
                _mm_add_epi64(_mm256_castsi256_si128(sum), _mm256_extracti64x2_epi64(sum, 1))
            };
            return _mm_cvtsi128_si64(sum2) + _mm_extract_epi64(sum2, 1);
        } else if constexpr (std::is_same_v<value_type, float>) {
            __m128 sum{_mm_mul_ps(
                _mm_maskz_load_ps(0b0111, reinterpret_cast<const float*>(&x)),
                _mm_maskz_load_ps(0b0111, reinterpret_cast<const float*>(&obj.x))
            )};
            sum = _mm_hadd_ps(sum, sum);
            return _mm_cvtss_f32(_mm_hadd_ps(sum, sum));
        } else if constexpr (std::is_same_v<value_type, double>) {
            const __m256d sum{_mm256_mul_pd(
                _mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&x)),
                _mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&obj.x))
            )};
            const __m128d sum2{
                _mm_add_pd(_mm256_castpd256_pd128(sum), _mm256_extractf64x2_pd(sum, 1))
            };
            return _mm_cvtsd_f64(_mm_hadd_pd(sum2, sum2));
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            const __m256 lhs{_mm256_maskz_load_ps(0b00111111, reinterpret_cast<const float*>(&x))};
            const __m256 rhs{
                _mm256_maskz_load_ps(0b00111111, reinterpret_cast<const float*>(&obj.x))
            };
            const __m256 sum{_mm256_addsub_ps(
                _mm256_mul_ps(lhs, _mm256_moveldup_ps(rhs)),
                _mm256_mul_ps(_mm256_permute_ps(lhs, 0xB1), _mm256_movehdup_ps(rhs))
            )};
            const __m128 sum2{
                _mm_add_ps(_mm256_castps256_ps128(sum), _mm256_extractf32x4_ps(sum, 1))
            };
            value_type result;
            _mm_storel_pi(
                reinterpret_cast<__m64*>(&result),
                _mm_add_ps(sum2, _mm_permute_ps(sum2, _MM_SHUFFLE(1, 0, 3, 2)))
            );
            return result;
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            const __m512d lhs{_mm512_maskz_load_pd(0b00111111, reinterpret_cast<const double*>(&x))
            };
            const __m512d rhs{
                _mm512_maskz_load_pd(0b00111111, reinterpret_cast<const double*>(&obj.x))
            };
            const __m512d sum{_mm512_add_pd(
                _mm512_mul_pd(lhs, _mm512_movedup_pd(rhs)),
                _mm512_xor_pd(
                    _mm512_mul_pd(_mm512_permute_pd(lhs, 0x55), _mm512_permute_pd(rhs, 0xFF)),
                    _mm512_set_pd(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0)
                )
            )};
            const __m256d sum2{
                _mm256_add_pd(_mm512_castpd512_pd256(sum), _mm512_extractf64x4_pd(sum, 1))
            };
            value_type result;
            _mm_store_pd(
                reinterpret_cast<double*>(&result),
                _mm_add_pd(_mm256_extractf64x2_pd(sum2, 0), _mm256_extractf64x2_pd(sum2, 1))
            );
            return result;
        } else {
            return x * obj.x + y * obj.y + z * obj.z;
        }
#else
        return x * obj.x + y * obj.y + z * obj.z;
#endif
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto cross(const Vector& obj) const noexcept -> Vector {
        Vector result;
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            const __m128i lhs{
                _mm_maskz_load_epi32(0b0111, reinterpret_cast<const std::int32_t*>(&x))
            };
            const __m128i rhs{
                _mm_maskz_load_epi32(0b0111, reinterpret_cast<const std::int32_t*>(&obj.x))
            };
            _mm_store_si128(
                reinterpret_cast<__m128i*>(&result.x),
                _mm_sub_epi32(
                    _mm_mullo_epi32(
                        _mm_shuffle_epi32(lhs, _MM_SHUFFLE(3, 0, 2, 1)),
                        _mm_shuffle_epi32(rhs, _MM_SHUFFLE(3, 1, 0, 2))
                    ),
                    _mm_mullo_epi32(
                        _mm_shuffle_epi32(lhs, _MM_SHUFFLE(3, 1, 0, 2)),
                        _mm_shuffle_epi32(rhs, _MM_SHUFFLE(3, 0, 2, 1))
                    )
                )
            );
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            const __m256i lhs{
                _mm256_maskz_load_epi64(0b0111, reinterpret_cast<const std::int64_t*>(&x))
            };
            const __m256i rhs{
                _mm256_maskz_load_epi64(0b0111, reinterpret_cast<const std::int64_t*>(&obj.x))
            };
            _mm256_store_si256(
                reinterpret_cast<__m256i*>(&result.x),
                _mm256_sub_epi64(
                    _mm256_mullo_epi64(
                        _mm256_permute4x64_epi64(lhs, _MM_SHUFFLE(3, 0, 2, 1)),
                        _mm256_permute4x64_epi64(rhs, _MM_SHUFFLE(3, 1, 0, 2))
                    ),
                    _mm256_mullo_epi64(
                        _mm256_permute4x64_epi64(lhs, _MM_SHUFFLE(3, 1, 0, 2)),
                        _mm256_permute4x64_epi64(rhs, _MM_SHUFFLE(3, 0, 2, 1))
                    )
                )
            );
        } else if constexpr (std::is_same_v<value_type, float>) {
            const __m128 lhs{_mm_maskz_load_ps(0b0111, reinterpret_cast<const float*>(&x))};
            const __m128 rhs{_mm_maskz_load_ps(0b0111, reinterpret_cast<const float*>(&obj.x))};
            _mm_store_ps(
                reinterpret_cast<float*>(&result.x),
                _mm_sub_ps(
                    _mm_mul_ps(
                        _mm_permute_ps(lhs, _MM_SHUFFLE(3, 0, 2, 1)),
                        _mm_permute_ps(rhs, _MM_SHUFFLE(3, 1, 0, 2))
                    ),
                    _mm_mul_ps(
                        _mm_permute_ps(lhs, _MM_SHUFFLE(3, 1, 0, 2)),
                        _mm_permute_ps(rhs, _MM_SHUFFLE(3, 0, 2, 1))
                    )
                )
            );
        } else if constexpr (std::is_same_v<value_type, double>) {
            const __m256d lhs{_mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&x))};
            const __m256d rhs{_mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&obj.x))
            };
            _mm256_store_pd(
                reinterpret_cast<double*>(&result.x),
                _mm256_sub_pd(
                    _mm256_mul_pd(
                        _mm256_permute4x64_pd(lhs, _MM_SHUFFLE(3, 0, 2, 1)),
                        _mm256_permute4x64_pd(rhs, _MM_SHUFFLE(3, 1, 0, 2))
                    ),
                    _mm256_mul_pd(
                        _mm256_permute4x64_pd(lhs, _MM_SHUFFLE(3, 1, 0, 2)),
                        _mm256_permute4x64_pd(rhs, _MM_SHUFFLE(3, 0, 2, 1))
                    )
                )
            );
        } else {
            result.x = y * obj.z - z * obj.y;
            result.y = z * obj.x - x * obj.z;
            result.z = x * obj.y - y * obj.x;
        }
#else
        result.x = y * obj.z - z * obj.y;
        result.y = z * obj.x - x * obj.z;
        result.z = x * obj.y - y * obj.x;
#endif
        return result;
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto crossProj(const Vector& obj) const noexcept -> value_type {
        return cross(obj).euclidean();
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto normalized() const noexcept -> Vector {
        return *this / euclidean();
    }
    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type> {
        return std::sqrt(euclideanSqr());
    }
    [[using gnu: pure, always_inline, leaf, hot]]
    constexpr auto euclideanSqr() const noexcept -> detail::DenseNormTraitT<value_type> {
#ifdef BOYLE_USE_AVX
        if constexpr (std::is_same_v<value_type, std::int32_t> ||
                      std::is_same_v<value_type, std::uint32_t>) {
            const __m128i lhs{
                _mm_maskz_load_epi32(0b0111, reinterpret_cast<const std::int32_t*>(&x))
            };
            __m128i sum{_mm_mullo_epi32(lhs, lhs)};
            sum = _mm_hadd_epi32(sum, sum);
            return _mm_cvtsi128_si32(_mm_hadd_epi32(sum, sum));
        } else if constexpr (std::is_same_v<value_type, std::int64_t> ||
                             std::is_same_v<value_type, std::uint64_t>) {
            const __m256i lhs{
                _mm256_maskz_load_epi64(0b0111, reinterpret_cast<const std::int64_t*>(&x))
            };
            const __m256i sum{_mm256_mullo_epi64(lhs, lhs)};
            const __m128i sum2{
                _mm_add_epi64(_mm256_castsi256_si128(sum), _mm256_extracti128_si256(sum, 1))
            };
            return _mm_cvtsi128_si64(sum2) + _mm_extract_epi64(sum2, 1);
        } else if constexpr (std::is_same_v<value_type, float>) {
            const __m128 lhs{_mm_maskz_load_ps(0b0111, reinterpret_cast<const float*>(&x))};
            __m128 sum{_mm_mul_ps(lhs, lhs)};
            sum = _mm_hadd_ps(sum, sum);
            return _mm_cvtss_f32(_mm_hadd_ps(sum, sum));
        } else if constexpr (std::is_same_v<value_type, double>) {
            const __m256d lhs{_mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&x))};
            const __m256d sum{_mm256_mul_pd(lhs, lhs)};
            const __m128d sum2{
                _mm_add_pd(_mm256_castpd256_pd128(sum), _mm256_extractf64x2_pd(sum, 1))
            };
            return _mm_cvtsd_f64(_mm_hadd_pd(sum2, sum2));
        } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
            const __m256 lhs{_mm256_maskz_load_ps(0b00111111, reinterpret_cast<const float*>(&x))};
            const __m256 sum{_mm256_mul_ps(lhs, lhs)};
            __m128 sum2{_mm_add_ps(_mm256_extractf32x4_ps(sum, 0), _mm256_extractf32x4_ps(sum, 1))};
            sum2 = _mm_hadd_ps(sum2, sum2);
            return _mm_cvtss_f32(_mm_hadd_ps(sum2, sum2));
        } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
            const __m512d lhs{_mm512_maskz_load_pd(0b00111111, reinterpret_cast<const double*>(&x))
            };
            const __m512d sum{_mm512_mul_pd(lhs, lhs)};
            const __m256d sum2{
                _mm256_add_pd(_mm512_extractf64x4_pd(sum, 0), _mm512_extractf64x4_pd(sum, 1))
            };
            const __m128d sum3{
                _mm_add_pd(_mm256_extractf64x2_pd(sum2, 0), _mm256_extractf64x2_pd(sum2, 1))
            };
            return _mm_cvtsd_f64(_mm_hadd_pd(sum3, sum3));
        } else {
            return std::norm(x) + std::norm(y) + std::norm(z);
        }
#else
        return std::norm(x) + std::norm(y) + std::norm(z);
#endif
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
        return dot(obj) < tol * tol;
    }

    [[using gnu: always_inline]]
    constexpr auto selfConjugated() noexcept -> Vector& {
#ifdef BOYLE_USE_AVX
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            if constexpr (std::is_same_v<value_type, std::complex<float>>) {
                _mm256_store_ps(
                    reinterpret_cast<float*>(&x),
                    _mm256_xor_ps(
                        _mm256_load_ps(reinterpret_cast<const float*>(&x)),
                        _mm256_set_ps(-0.0F, 0.0F, -0.0F, 0.0F, -0.0F, 0.0F, -0.0F, 0.0F)
                    )
                );
            } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
                _mm512_store_pd(
                    reinterpret_cast<double*>(&x),
                    _mm512_xor_pd(
                        _mm512_load_pd(reinterpret_cast<const double*>(&x)),
                        _mm512_set_pd(-0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0)
                    )
                );
            } else {
                x = std::conj(x);
                y = std::conj(y);
                z = std::conj(z);
            }
        }
#else
        if constexpr (::boyle::math::isComplexArithmeticV<value_type>) {
            x = std::conj(x);
            y = std::conj(y);
            z = std::conj(z);
        }
#endif
        return *this;
    }
    [[using gnu: pure, always_inline]]
    constexpr auto conjugated() const& noexcept -> Vector {
        Vector result{*this};
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

    value_type x, y, z;
};

template <ScalarArithmetic T>
[[using gnu: pure, always_inline, hot]]
inline constexpr auto operator*(const ScalarArithmetic auto& fac, const Vector<T, 3>& obj) noexcept
    -> Vector<T, 3> {
    return obj * fac;
}

template <ScalarArithmetic T>
[[using gnu: always_inline, hot]]
inline constexpr auto operator*(const ScalarArithmetic auto& fac, Vector<T, 3>&& obj) noexcept
    -> Vector<T, 3>&& {
    return std::move(obj) * fac;
}

template <typename Char, ScalarArithmetic T>
[[using gnu: always_inline]]
inline auto operator<<(std::basic_ostream<Char>& os, const Vector<T, 3>& obj) noexcept
    -> std::basic_ostream<Char>& {
    os << "(x: " << obj.x << ", y: " << obj.y << ", z: " << obj.z << ")";
    return os;
}

template <ScalarArithmetic T>
using Vec3 = Vector<T, 3>;

using Vec3s = Vec3<float>;
using Vec3d = Vec3<double>;
using Vec3c = Vec3<std::complex<float>>;
using Vec3z = Vec3<std::complex<double>>;

template <typename T>
struct isVec3Arithmetic : std::false_type {};

template <ScalarArithmetic T>
struct isVec3Arithmetic<Vec3<T>> : std::true_type {};

template <typename T>
inline constexpr bool isVec3ArithmeticV = isVec3Arithmetic<T>::value;

template <typename T>
concept Vec3Arithmetic = isVec3ArithmeticV<T>;

template <Vec3Arithmetic T>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline auto squeeze(
    std::ranges::forward_range auto&& xs, std::ranges::forward_range auto&& ys,
    std::ranges::forward_range auto&& zs
) noexcept -> std::vector<T>
    requires std::same_as<std::ranges::range_value_t<decltype(xs)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(ys)>, typename T::value_type> &&
             std::same_as<std::ranges::range_value_t<decltype(zs)>, typename T::value_type>
{
    if (xs.size() != ys.size() || ys.size() != zs.size()) [[unlikely]] {
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

template <std::floating_point T, typename Char>
struct formatter<::boyle::math::Vec3<T>, Char> {
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
    auto format(const ::boyle::math::Vec3<T>& obj, FormatContext& ctx) const
        -> decltype(ctx.out()) {
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

template <std::floating_point T>
[[using gnu: pure, always_inline]]
inline constexpr auto hypot(const ::boyle::math::Vec3<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return std::hypot(obj.x, obj.y, obj.z);
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto abs(const ::boyle::math::Vec3<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return obj.euclidean();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto norm(const ::boyle::math::Vec3<T>& obj) noexcept
    -> ::boyle::math::detail::DenseNormTraitT<T> {
    return obj.euclideanSqr();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: pure, always_inline]]
inline constexpr auto conj(const ::boyle::math::Vec3<T>& vector) noexcept
    -> ::boyle::math::Vec3<T> {
    return vector.conjugated();
}

template <::boyle::math::ScalarArithmetic T>
[[using gnu: always_inline]]
inline constexpr auto conj(::boyle::math::Vec3<T>&& vector) noexcept -> ::boyle::math::Vec3<T>&& {
    return std::move(vector).conjugated();
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)

namespace fmt {

template <::boyle::math::ScalarArithmetic T, typename Char>
struct formatter<::boyle::math::Vec3<T>, Char> {
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
    auto format(const ::boyle::math::Vec3<T>& obj, FormatContext& ctx) const
        -> decltype(ctx.out()) {
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

template <::boyle::math::ScalarArithmetic T>
[[using gnu: always_inline]]
inline constexpr auto serialize(
    auto& archive, ::boyle::math::Vec3<T>& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.x;
    archive & obj.y;
    archive & obj.z;
    return;
}

} // namespace boost::serialization
