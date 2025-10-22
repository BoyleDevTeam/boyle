/**
 * @file vector_view.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-04-08
 *
 * @copyright Copyright (c) 2025 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <type_traits>

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/dense_traits.hpp"
#include "boyle/math/dense/detail/dense_norm_trait.hpp"
#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <ScalarArithmetic Scalar>
class VectorView final {
  public:
    using value_type = typename DenseTraits<VectorView>::value_type;
    using reference = typename DenseTraits<VectorView>::reference;
    using const_reference = typename DenseTraits<VectorView>::const_reference;
    using pointer = typename DenseTraits<VectorView>::pointer;
    using const_pointer = typename DenseTraits<VectorView>::const_pointer;
    using size_type = typename DenseTraits<VectorView>::size_type;
    using difference_type = typename DenseTraits<VectorView>::difference_type;

    constexpr VectorView() noexcept = delete;
    constexpr VectorView(const VectorView& other) noexcept = default;
    constexpr auto operator=(const VectorView& other) noexcept -> VectorView& = default;
    constexpr VectorView(VectorView&& other) noexcept = default;
    constexpr auto operator=(VectorView&& other) noexcept -> VectorView& = default;
    constexpr ~VectorView() noexcept = default;

    [[using gnu: always_inline]]
    explicit VectorView(pointer data, size_type size, size_type stride = 1) noexcept
        : m_data{data}, m_size{size}, m_stride{stride} {}

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto size() const noexcept -> size_type {
        return m_size;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto stride() const noexcept -> size_type {
        return m_stride;
    }

    [[using gnu: pure, always_inline]]
    auto data() noexcept -> pointer {
        return m_data;
    }
    [[using gnu: pure, always_inline]]
    auto data() const noexcept -> const_pointer {
        return m_data;
    }

    [[using gnu: always_inline, hot]]
    constexpr auto fill(const_reference value) noexcept -> void {
        for (size_type i{0}; i < m_size; ++i) {
            operator[](i) = value;
        }
        return;
    }

    [[nodiscard]]
    auto empty() const noexcept -> bool {
        return m_size == 0;
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type i) noexcept -> reference {
        return m_data[i * m_stride];
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator[](size_type i) const noexcept -> const_reference {
        return m_data[i * m_stride];
    }

    [[using gnu: pure, always_inline, hot]] [[nodiscard]]
    constexpr auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Vector index out of range.");
        }
#endif
        return operator[](i);
    }

    [[using gnu: always_inline, hot]]
    constexpr auto updateCoeff(size_type i, const_reference value) const noexcept(
        !BOYLE_CHECK_PARAMS
    ) -> void {
#if BOYLE_CHECK_PARAMS == 1
        if (i >= size()) {
            throw std::out_of_range("Vector index out of range.");
        }
#endif
        operator[](i) = value;
        return;
    }

    [[using gnu: pure, always_inline]]
    auto view(std::pair<size_type, size_type> range) noexcept -> VectorView<value_type> {
        range.second = std::min(range.second, size());
        return VectorView<value_type>{
            data() + range.first * stride(), range.second - range.first, stride()
        };
    }

    [[using gnu: pure, always_inline]]
    auto view(std::pair<size_type, size_type> range) const noexcept
        -> VectorView<const value_type> {
        range.second = std::min(range.second, size());
        return VectorView<value_type>{
            data() + range.first * stride(), range.second - range.first, stride()
        };
    }

    [[using gnu: pure, always_inline]]
    auto view() noexcept -> VectorView<value_type> {
        return view({0, std::numeric_limits<value_type>::max()});
    }

    [[using gnu: pure, always_inline]]
    auto view() const noexcept -> VectorView<const value_type> {
        return view({0, std::numeric_limits<value_type>::max()});
    }

    [[using gnu: pure, always_inline, hot]]
    constexpr auto operator==(const VectorView& other) const noexcept -> bool {
        if (size() != other.size()) [[unlikely]] {
            return false;
        }
        const size_type n{size()};
        for (size_type i{0}; i < n; ++i) {
            if (m_data[i] != other.m_data[i]) {
                return false;
            }
        }
        return true;
    }

    [[using gnu: pure, always_inline, hot]]
    auto identicalTo(VectorView<value_type> obj, double tol = 1E-8) const noexcept -> bool {
        if (size() != obj.size()) {
            return false;
        }
        size_type n{size()};
        detail::DenseNormTraitT<value_type> test(0.0);
        for (size_type i{0}; i < n; ++i) {
            test = std::max(
                test, std::abs(operator[](i) - obj[i]) /
                          std::max(std::abs(operator[](i)), static_cast<decltype(test)>(1.0))
            );
        }
        return test < tol;
    }

  private:
    pointer m_data;
    size_type m_size, m_stride;
};

template <typename Char, ScalarArithmetic Scalar>
inline auto operator<<(std::basic_ostream<Char>& os, VectorView<Scalar> vector_view) noexcept
    -> std::basic_ostream<Char>& {
    using value_type = std::remove_const_t<Scalar>;
    using size_type = typename VectorView<Scalar>::size_type;
    constexpr size_type kWidth{isComplexArithmeticV<value_type> ? 32 : 16};
    const size_type size{vector_view.size()};
    os << std::fixed;
    for (size_type i{0}; i < size; ++i) {
        os << std::setw(kWidth) << vector_view[i] << '\n';
    }
    return os;
}

} // namespace boyle::math
