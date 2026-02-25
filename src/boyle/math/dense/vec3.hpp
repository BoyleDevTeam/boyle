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

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <format>
#include <ostream>

#if BOYLE_USE_SIMD == 1
#include <immintrin.h>
#endif

#include "boost/serialization/complex.hpp"

#include "boyle/common/utils/in_in_in_out_result.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic T>
class alignas(4 * sizeof(T)) Vector<T, 3> final {
  public:
    using value_type = typename DenseTraits<Vector>::value_type;
    using reference = typename DenseTraits<Vector>::reference;
    using const_reference = typename DenseTraits<Vector>::const_reference;
    using pointer = typename DenseTraits<Vector>::pointer;
    using const_pointer = typename DenseTraits<Vector>::const_pointer;
    using size_type = typename DenseTraits<Vector>::size_type;
    using difference_type = typename DenseTraits<Vector>::difference_type;
    using allocator_type = typename DenseTraits<Vector>::allocator_type;

    static constexpr size_type kSize = 2;

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

    [[using gnu: always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return allocator_type{};
    }

    [[using gnu: always_inline]]
    explicit Vector([[maybe_unused]] const allocator_type& alloc = {}) noexcept {}

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

    [[using gnu: const, always_inline, leaf]]
    static constexpr auto size() noexcept -> size_type {
        return kSize;
    }
    [[using gnu: const, always_inline, leaf]]
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

    [[using gnu: pure, always_inline, hot]]
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
#if BOYLE_USE_SIMD == 1
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
            const __m512d lhs{
                _mm512_maskz_load_pd(0b00111111, reinterpret_cast<const double*>(&x))
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

    template <size_type OuterSize, MatrixOrder ObjOrder>
    [[using gnu: pure, always_inline, hot]]
    constexpr auto dot(const Matrix<value_type, kSize, OuterSize, ObjOrder>& obj) const noexcept
        -> Vector {
        Vector<value_type, OuterSize> result(OuterSize, 0.0);
        for (size_type j{0}; j < OuterSize; ++j) {
            for (size_type i{0}; i < kSize; ++i) {
                result[j] += operator[](i) * obj[i, j];
            }
        }
        return result;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto cross(const Vector& obj) const noexcept -> Vector {
        Vector result;
#if BOYLE_USE_SIMD == 1
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
            const __m256d rhs{
                _mm256_maskz_load_pd(0b0111, reinterpret_cast<const double*>(&obj.x))
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
#if BOYLE_USE_SIMD == 1
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
            const __m512d lhs{
                _mm512_maskz_load_pd(0b00111111, reinterpret_cast<const double*>(&x))
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
#if BOYLE_USE_SIMD == 1
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
    obj *= fac;
    return std::move(obj);
}

template <typename Char, ScalarArithmetic T>
[[using gnu: always_inline]]
inline auto operator<<(std::basic_ostream<Char>& os, const Vector<T, 3>& obj)
    -> std::basic_ostream<Char>& {
    os << std::format("{}", obj);
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

template <
    std::input_iterator InputIt1, std::sentinel_for<InputIt1> SentinelIt1,
    std::input_iterator InputIt2, std::sentinel_for<InputIt1> SentinelIt2,
    std::input_iterator InputIt3, std::sentinel_for<InputIt1> SentinelIt3,
    std::weakly_incrementable OutputIt>
    requires std::floating_point<typename std::iterator_traits<InputIt1>::value_type> &&
             std::floating_point<typename std::iterator_traits<InputIt2>::value_type> &&
             std::floating_point<typename std::iterator_traits<InputIt3>::value_type> &&
             Vec3Arithmetic<typename std::iterator_traits<OutputIt>::value_type>
[[using gnu: always_inline, hot]]
inline auto squeeze(
    InputIt1 first1, SentinelIt1 last1, InputIt2 first2, SentinelIt2 last2, InputIt3 first3,
    SentinelIt3 last3, OutputIt result
) noexcept -> ::boyle::common::InInInOutResult<InputIt1, InputIt2, InputIt3, OutputIt> {
    for (; first1 != last1 && first2 != last2 && first3 != last3;
         ++first1, ++first2, ++first3, ++result) {
        *result = typename std::iterator_traits<OutputIt>::value_type{*first1, *first2, *first3};
    }
    return {first1, first2, first3, result};
}

template <
    std::ranges::input_range InputRange1, std::ranges::input_range InputRange2,
    std::ranges::input_range InputRange3, std::weakly_incrementable OutputIt>
    requires std::floating_point<typename std::ranges::range_value_t<InputRange1>> &&
             std::floating_point<typename std::ranges::range_value_t<InputRange2>> &&
             std::floating_point<typename std::ranges::range_value_t<InputRange3>> &&
             Vec3Arithmetic<typename std::iterator_traits<OutputIt>::value_type>
[[using gnu: always_inline, hot]]
inline auto squeeze(InputRange1&& xs, InputRange2&& ys, InputRange3&& zs, OutputIt result) noexcept
    -> ::boyle::common::InInInOutResult<
        std::ranges::iterator_t<InputRange1>, std::ranges::iterator_t<InputRange2>,
        std::ranges::iterator_t<InputRange3>, OutputIt> {
    return squeeze(
        std::ranges::begin(xs), std::ranges::end(xs), std::ranges::begin(ys), std::ranges::end(ys),
        std::ranges::begin(zs), std::ranges::end(zs), result
    );
}

} // namespace boyle::math

// NOLINTBEGIN(cert-dcl58-cpp)

namespace std {

template <std::floating_point T, typename Char>
struct formatter<::boyle::math::Vec3<T>, Char> {
    [[using gnu: ]]
    constexpr auto parse(std::basic_format_parse_context<Char>& ctx) -> decltype(ctx.begin()) {
        auto it{ctx.begin()};
        const auto end{ctx.end()};
        if (it != end && *it >= static_cast<Char>('0') && *it <= static_cast<Char>('9')) {
            width = 0;
            for (; it != end && *it >= static_cast<Char>('0') && *it <= static_cast<Char>('9');
                 ++it) {
                width = width * 10 + (*it - static_cast<Char>('0'));
            }
        }
        if (it != end && *it == static_cast<Char>('.')) {
            ++it;
            precision = 0;
            for (; it != end && *it >= static_cast<Char>('0') && *it <= static_cast<Char>('9');
                 ++it) {
                precision = precision * 10 + (*it - static_cast<Char>('0'));
            }
        }
        return it;
    }

    template <typename FormatContext>
    [[using gnu: ]]
    auto format(const ::boyle::math::Vec3<T>& obj, FormatContext& ctx) const
        -> decltype(ctx.out()) {
        auto out = ctx.out();
        if (width > 0 && precision >= 0) {
            return std::format_to(
                out, "(x: {:>{}.{}f}, y: {:>{}.{}f}, z: {:>{}.{}f})", obj.x, width, precision,
                obj.y, width, precision, obj.z, width, precision
            );
        }
        if (width > 0) {
            return std::format_to(
                out, "(x: {:>{}f}, y: {:>{}f}, z: {:>{}f})", obj.x, width, obj.y, width, obj.z,
                width
            );
        }
        if (precision >= 0) {
            return std::format_to(
                out, "(x: {:.{}f}, y: {:.{}f}, z: {:.{}f})", obj.x, precision, obj.y, precision,
                obj.z, precision
            );
        }
        return std::format_to(out, "(x: {}, y: {}, z: {})", obj.x, obj.y, obj.z);
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
    vector.selfConjugated();
    return std::move(vector);
}

} // namespace std

// NOLINTEND(cert-dcl58-cpp)

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
