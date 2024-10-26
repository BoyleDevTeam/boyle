/**
 * @file chebyshev.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-02-06
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>

#include "boost/serialization/access.hpp"
#include "boost/serialization/array.hpp"

namespace boyle::math {

template <std::floating_point Scalar = double, std::size_t Size = 10>
class [[nodiscard]] Chebyshev final {
    friend class boost::serialization::access;

  public:
    Chebyshev() noexcept = default;
    Chebyshev(const Chebyshev&) noexcept = default;
    auto operator=(const Chebyshev&) noexcept -> Chebyshev& = default;
    Chebyshev(Chebyshev&&) noexcept = default;
    auto operator=(Chebyshev&&) noexcept -> Chebyshev& = default;
    ~Chebyshev() noexcept = default;

    [[using gnu: ]]
    explicit Chebyshev(
        std::function<auto(Scalar)->Scalar> func, Scalar lower_bound, Scalar upper_bound,
        Scalar abs_tol = 1E-8
    ) noexcept
        : m_lower_bound{lower_bound}, m_upper_bound{upper_bound} {
        constexpr Scalar kFactor{2.0 / Size};
        const Scalar middle{(m_upper_bound + m_lower_bound) * 0.5};
        const Scalar radius{(m_upper_bound - m_lower_bound) * 0.5};

        std::array<Scalar, Size> func_knots;
        for (std::size_t i{0}; i < Size; ++i) {
            func_knots[i] = func(middle + radius * std::cos(M_PI * (i + 0.5) / Size));
        }

        for (std::size_t i{0}; i < Size; ++i) {
            m_coeffs[i] = 0.0;
            for (std::size_t j{0}; j < Size; ++j) {
                m_coeffs[i] += func_knots[j] * std::cos(M_PI * i * (j + 0.5) / Size);
            }
            m_coeffs[i] *= kFactor;
        }

        truncate(abs_tol);
    }

    template <std::floating_point OtherScalar, std::size_t OtherSize>
    [[using gnu: always_inline, leaf]]
    explicit Chebyshev(const Chebyshev<OtherScalar, OtherSize>& other) noexcept
        requires(Size >= OtherSize)
    {
        m_lower_bound = other.m_lower_bound;
        m_upper_bound = other.m_upper_bound;
        m_num_trunc = other.m_num_trunc;
        for (std::size_t i{0}; i < OtherSize; ++i) {
            m_coeffs[i] = other.m_coeffs[i];
        }
    }

    [[using gnu: always_inline]]
    auto truncate(Scalar abs_tol) noexcept -> void {
        m_num_trunc = Size;
        while (m_num_trunc > 1 && std::abs(m_coeffs[m_num_trunc - 1]) < abs_tol) {
            --m_num_trunc;
        }
        return;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    constexpr auto size() const noexcept -> std::size_t {
        return Size;
    }

    [[using gnu: pure, always_inline, leaf]] [[nodiscard]]
    auto num_trunc() const noexcept -> std::size_t {
        return m_num_trunc;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto coeff(std::size_t i) const noexcept -> Scalar {
        return m_coeffs[i];
    }

    [[using gnu: pure, always_inline, leaf]]
    auto coeffs() const noexcept -> const std::array<Scalar, Size>& {
        return m_coeffs;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto lower_bound() const noexcept -> Scalar {
        return m_lower_bound;
    }

    [[using gnu: pure, always_inline, leaf]]
    auto upper_bound() const noexcept -> Scalar {
        return m_upper_bound;
    }

    [[using gnu: pure, always_inline]]
    auto eval(Scalar t) const noexcept -> Scalar {
        return eval(t, m_num_trunc);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto eval(Scalar t, std::size_t num_trunc) const noexcept -> Scalar {
        if ((t - m_lower_bound) * (t - m_upper_bound) > 0.0) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        const Scalar y{(t * 2.0 - m_lower_bound - m_upper_bound) / (m_upper_bound - m_lower_bound)};
        const Scalar y2{y * 2.0};
        Scalar d{0.0}, dd{0.0};
        for (int i = num_trunc - 1; i > 0; --i) {
            const Scalar sv{d};
            d = y2 * d - dd + m_coeffs[i];
            dd = sv;
        }

        return y * d - dd + 0.5 * m_coeffs[0];
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(Scalar t) const noexcept -> Scalar {
        return eval(t);
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(Scalar t, std::size_t num_trunc) const noexcept -> Scalar {
        return eval(t, num_trunc);
    }

    [[using gnu: pure]]
    auto derivative() const noexcept -> Chebyshev {
        const Scalar factor{2.0 / (m_upper_bound - m_lower_bound)};
        Chebyshev derivative;
        derivative.m_lower_bound = m_lower_bound;
        derivative.m_upper_bound = m_upper_bound;
        derivative.m_num_trunc = m_num_trunc;

        derivative.m_coeffs[m_num_trunc - 1] = 0.0;
        derivative.m_coeffs[m_num_trunc - 2] = m_coeffs[m_num_trunc - 1] * (m_num_trunc - 1) * 2;
        for (int i = m_num_trunc - 2; i > 0; --i) {
            derivative.m_coeffs[i - 1] = derivative.m_coeffs[i + 1] + m_coeffs[i] * i * 2;
        }
        for (Scalar& coeff : derivative.m_coeffs) {
            coeff *= factor;
        }

        return derivative;
    }

    [[using gnu: pure]]
    auto integral() const noexcept -> Chebyshev {
        const Scalar factor{(m_upper_bound - m_lower_bound) * 0.25};
        Chebyshev integral;
        integral.m_lower_bound = m_lower_bound;
        integral.m_upper_bound = m_upper_bound;
        integral.m_num_trunc = m_num_trunc;

        Scalar sum{0.0}, sign{1.0};
        for (std::size_t i{1}; i < m_num_trunc - 1; ++i) {
            integral.m_coeffs[i] = factor * (m_coeffs[i - 1] - m_coeffs[i + 1]) / i;
            sum += sign * integral.m_coeffs[i];
            sign = -sign;
        }
        integral.m_coeffs[m_num_trunc - 1] = factor * m_coeffs[m_num_trunc - 2] / (m_num_trunc - 1);
        sum += sign * integral.m_coeffs[m_num_trunc - 1];
        integral.m_coeffs[0] = sum * 2.0;

        return integral;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_coeffs;
        archive & m_lower_bound;
        archive & m_upper_bound;
        archive & m_num_trunc;
        return;
    }

    std::array<Scalar, Size> m_coeffs{};
    Scalar m_lower_bound{0.0}, m_upper_bound{1.0};
    std::size_t m_num_trunc{0};
};

} // namespace boyle::math
