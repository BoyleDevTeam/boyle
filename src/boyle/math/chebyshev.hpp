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

#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"

namespace boyle::math {

template <std::floating_point Scalar = double, std::integral Index = int>
class [[nodiscard]] Chebyshev final {
    friend class boost::serialization::access;

  public:
    [[using gnu: flatten, leaf]]
    explicit Chebyshev(
        std::function<auto(Scalar)->Scalar> func, Scalar lower_bound, Scalar upper_bound,
        Index num_total = 50
    ) noexcept
        : m_lower_bound{lower_bound}, m_upper_bound{upper_bound}, m_num_trunc{num_total} {
        const Scalar factor{2.0 / num_total};
        const Scalar upplo{(m_upper_bound + m_lower_bound) * 0.5};
        const Scalar upmlo{(m_upper_bound - m_lower_bound) * 0.5};

        std::vector<Scalar> func_knots;
        func_knots.reserve(num_total);

        for (Index i{0}; i < num_total; ++i) {
            func_knots.push_back(func(upplo + upmlo * std::cos(M_PI * (i + 0.5) / num_total)));
        }

        m_coeffs.reserve(num_total);
        for (Index i{0}; i < num_total; ++i) {
            Scalar sum{0.0};
            for (Index j{0}; j < num_total; ++j) {
                sum += func_knots[j] * std::cos(M_PI * i * (j + 0.5) / num_total);
            }
            m_coeffs.push_back(factor * sum);
        }
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(Chebyshev);

    ~Chebyshev() noexcept = default;

    [[using gnu: always_inline]]
    auto truncate(Scalar abs_tol) noexcept -> void {
        while (m_num_trunc > 1 && std::abs(m_coeffs[m_num_trunc - 1]) < abs_tol) {
            --m_num_trunc;
        }
        return;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto num_total() const noexcept -> std::size_t {
        return m_coeffs.size();
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto num_trunc() const noexcept -> std::size_t {
        return static_cast<std::size_t>(m_num_trunc);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto coeffs() const noexcept -> const std::vector<Scalar>& {
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
    auto operator()(Scalar t) const noexcept -> Scalar {
        return operator()(t, m_num_trunc);
    }

    [[using gnu: pure, always_inline, leaf]]
    auto operator()(Scalar t, Index num_trunc) const noexcept -> Scalar {
        if ((t - m_lower_bound) * (t - m_upper_bound) > 0.0) {
            BOYLE_LOG_ERROR(
                "Invalid argument detected! the variable t has to be in range of [m_lower_bound, "
                "m_upper_bound]: t = {0:f}, m_lower_bound = {1:f}, m_upper_bound = {2:f}",
                t, m_lower_bound, m_upper_bound
            );
            return std::numeric_limits<Scalar>::quiet_NaN();
        }

        const Scalar y{(t * 2.0 - m_lower_bound - m_upper_bound) / (m_upper_bound - m_lower_bound)};
        const Scalar y2{y * 2.0};
        Scalar d{0.0}, dd{0.0};
        for (int i{num_trunc - 1}; i > 0; --i) {
            const Scalar sv{d};
            d = y2 * d - dd + m_coeffs[i];
            dd = sv;
        }
        return y * d - dd + 0.5 * m_coeffs[0];
    }

    [[using gnu: pure, always_inline, leaf]]
    auto derivative() const noexcept -> Chebyshev {
        const Scalar factor{2.0 / (m_upper_bound - m_lower_bound)};
        const std::size_t num_total{m_coeffs.size()};
        std::vector<Scalar> dcoeffs(num_total);

        dcoeffs[num_total - 1] = 0.0;
        dcoeffs[num_total - 2] = m_coeffs[num_total - 1] * (num_total - 1) * 2;
        for (int i = num_total - 2; i > 0; --i) {
            dcoeffs[i - 1] = dcoeffs[i + 1] + m_coeffs[i] * i * 2;
        }

        for (Scalar& dcoeff : dcoeffs) {
            dcoeff *= factor;
        }
        return Chebyshev{std::move(dcoeffs), m_lower_bound, m_upper_bound};
    }

    [[using gnu: pure, always_inline, leaf]]
    auto integral() const noexcept -> Chebyshev {
        const Scalar factor{(m_upper_bound - m_lower_bound) * 0.25};
        const std::size_t num_total{m_coeffs.size()};
        std::vector<Scalar> integral_coeffs(num_total);

        Scalar sum{0.0}, sign{1.0};
        for (std::size_t i{1}; i < num_total - 1; ++i) {
            integral_coeffs[i] = factor * (m_coeffs[i - 1] - m_coeffs[i + 1]) / i;
            sum += sign * integral_coeffs[i];
            sign = -sign;
        }
        integral_coeffs[num_total - 1] = factor * m_coeffs[num_total - 2] / (num_total - 1);
        sum += sign * integral_coeffs[num_total - 1];
        integral_coeffs[0] = sum * 2.0;
        return Chebyshev(std::move(integral_coeffs), m_lower_bound, m_upper_bound);
    }

  private:
    [[using gnu: always_inline, leaf]]
    explicit Chebyshev(std::vector<Scalar> coeffs, Scalar lower_bound, Scalar upper_bound) noexcept
        : m_coeffs{std::move(coeffs)}, m_lower_bound{lower_bound}, m_upper_bound{upper_bound},
          m_num_trunc{static_cast<Index>(m_coeffs.size())} {}

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_coeffs;
        archive & m_lower_bound;
        archive & m_upper_bound;
        archive & m_num_trunc;
        return;
    }

    std::vector<Scalar> m_coeffs{};
    Scalar m_lower_bound{0.0}, m_upper_bound{1.0};
    Index m_num_trunc{0};
};

} // namespace boyle::math
