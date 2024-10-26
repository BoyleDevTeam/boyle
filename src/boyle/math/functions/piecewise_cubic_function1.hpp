/**
 * @file piecewise_cubic_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "fmt/format.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/cubic_interpolation.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"
#include "boyle/math/vec3.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] PiecewiseCubicFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr param_type kDuplicateCriterion{1E-8};

    struct BoundaryMode final {
        unsigned int order;
        value_type derivative;
    };

    PiecewiseCubicFunction1() noexcept = default;
    PiecewiseCubicFunction1(const PiecewiseCubicFunction1& other) noexcept = default;
    auto operator=(const PiecewiseCubicFunction1& other
    ) noexcept -> PiecewiseCubicFunction1& = default;
    PiecewiseCubicFunction1(PiecewiseCubicFunction1&& other) noexcept = default;
    auto operator=(PiecewiseCubicFunction1&& other) noexcept -> PiecewiseCubicFunction1& = default;
    ~PiecewiseCubicFunction1() noexcept = default;

    [[using gnu: always_inline]]
    explicit PiecewiseCubicFunction1(std::vector<param_type> ts, std::vector<value_type> ys)
        : PiecewiseCubicFunction1(
              std::move(ts), std::move(ys), BoundaryMode{2, 0.0}, BoundaryMode{2, 0.0}
          ) {}

    [[using gnu: ]]
    explicit PiecewiseCubicFunction1(
        std::vector<param_type> ts, std::vector<value_type> ys, BoundaryMode b0, BoundaryMode bf
    ) noexcept(!BOYLE_CHECK_PARAMS)
        : m_ts{std::move(ts)}, m_ys{std::move(ys)} {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ts.size() < 2 || m_ys.size() < 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                m_ts.size(), m_ys.size()
            ));
        }
        if (m_ts.size() != m_ys.size()) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                m_ts.size(), m_ys.size()
            ));
        }
        if (!std::is_sorted(m_ts.cbegin(), m_ts.cend())) [[unlikely]] {
            throw std::invalid_argument(
                fmt::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        }
        if (hasDuplicates(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, kDuplicateCriterion))
            [[unlikely]] {
            throw std::invalid_argument(
                fmt::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }

        if (b0.order == 0 || b0.order > 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "b0.order = {0:d}.",
                b0.order
            ));
        }
        if (bf.order == 0 || bf.order > 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "bf.order = {0:d}.",
                bf.order
            ));
        }
#endif
        const std::size_t size{m_ts.size()};
        std::vector<param_type> hs(size - 1);
        std::vector<value_type> ds(size - 1);
        std::vector<param_type> a_diag(size);
        std::vector<value_type> b(size);

        hs[0] = m_ts[1] - m_ts[0];
        ds[0] = (m_ys[1] - m_ys[0]) / hs[0];
        for (std::size_t i{1}; i < size - 1; i++) {
            hs[i] = m_ts[i + 1] - m_ts[i];
            ds[i] = (m_ys[i + 1] - m_ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }

        std::vector<param_type> a_low{hs};
        std::vector<param_type> a_up{hs};

        if (b0.order == 2) {
            a_diag[0] = 1.0;
            a_up[0] = 0.0;
            b[0] = b0.derivative;
        } else {
            a_diag[0] = 1.0;
            a_up[0] = 0.5;
            b[0] = (ds[0] - b0.derivative) * 3.0 / hs[0];
        }
        if (bf.order == 2) {
            a_diag[size - 1] = 1.0;
            a_low[size - 2] = 0.0;
            b[size - 1] = bf.derivative;
        } else {
            a_diag[size - 1] = 1.0;
            a_low[size - 2] = 0.5;
            b[size - 1] = (bf.derivative - ds[size - 2]) * 3.0 / hs[size - 2];
        }

        const TridiagonalMatrix A{std::move(a_low), std::move(a_diag), std::move(a_up)};

        m_ddys = A.luDcmp(b);
    }

    [[using gnu: ]]
    explicit PiecewiseCubicFunction1(
        [[maybe_unused]] periodic_tag tag, std::vector<param_type> ts, std::vector<value_type> ys
    ) noexcept(!BOYLE_CHECK_PARAMS)
        : m_ts{std::move(ts)}, m_ys{std::move(ys)} {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ts.size() < 2 || m_ys.size() < 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                m_ts.size(), m_ys.size()
            ));
        }
        if (m_ts.size() != m_ys.size()) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                m_ts.size(), m_ys.size()
            ));
        }
        if (!std::is_sorted(m_ts.cbegin(), m_ts.cend())) [[unlikely]] {
            throw std::invalid_argument(
                fmt::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        }
        if (hasDuplicates(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, kDuplicateCriterion))
            [[unlikely]] {
            throw std::invalid_argument(
                fmt::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }
        if (std::abs(m_ys.front() - m_ys.back()) > kDuplicateCriterion) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! When choosing periodic boundary condition, it "
                "requires ys.front() == m_ys.back() while m_ys.front() = {0}, m_ys.back() = {1} "
                "here.",
                m_ys.front(), m_ys.back()
            ));
        }
#endif
        const std::size_t size{m_ts.size() - 1};
        std::vector<param_type> hs(size);
        std::vector<value_type> ds(size);
        std::vector<param_type> a_diag(size);
        std::vector<value_type> b(size);

        hs[0] = m_ts[1] - m_ts[0];
        ds[0] = (m_ys[1] - m_ys[0]) / hs[0];
        for (std::size_t i{1}; i < size; i++) {
            hs[i] = m_ts[i + 1] - m_ts[i];
            ds[i] = (m_ys[i + 1] - m_ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }
        a_diag[0] = (hs[0] + hs[size - 1]) * 2.0;
        b[0] = (ds[0] - ds[size - 1]) * 6.0;

        std::vector<param_type> a_low{hs.cbegin(), hs.cend() - 1};
        std::vector<param_type> a_up{hs.cbegin(), hs.cend() - 1};
        const param_type a_bottom{hs[size - 1]};
        const param_type a_top{hs[size - 1]};

        const PeriodicTridiagonalMatrix A{
            a_bottom, std::move(a_low), std::move(a_diag), std::move(a_up), a_top
        };

        m_ddys = A.luDcmp(b);
        m_ddys.push_back(m_ddys[0]);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto eval(param_type t) const noexcept -> value_type {
        constexpr std::array<param_type, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            const param_type h{m_ts[1] - m_ts[0]};
            const param_type ratio{(t - m_ts[0]) / h};
            return lerp(m_ys[0], m_ys[1], ratio) +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * (t - m_ts[0]) * h;
        }
        if (pos == m_ts.size()) {
            const param_type h{m_ts[pos - 2] - m_ts[pos - 1]};
            const param_type ratio{(t - m_ts[pos - 1]) / h};
            return lerp(m_ys[pos - 1], m_ys[pos - 2], ratio) +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) *
                       (t - m_ts[pos - 1]) * h;
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return cuberp(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative(param_type t) const noexcept -> value_type {
        constexpr std::array<param_type, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            const param_type h{m_ts[1] - m_ts[0]};
            return (m_ys[1] - m_ys[0]) / h +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * h;
        }
        if (pos == m_ts.size()) {
            const param_type h{m_ts[pos - 2] - m_ts[pos - 1]};
            return (m_ys[pos - 2] - m_ys[pos - 1]) / h +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) * h;
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio{(t - m_ts[pos - 1]) / h};
        return cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    auto derivative(param_type t, unsigned int order) const
        noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (order < 1 || order > 3) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid argument error! The PiecewiseLinearFunction only has 1, 2, 3 order "
                "derivatives: order = {0:d}.",
                order
            ));
        }
#endif
        value_type result;
        switch (order) {
        case 1:
            result = derivative(t);
            break;
        case 2:
            result = derivative2(t);
            break;
        case 3:
            result = derivative3(t);
            break;
        default:
            std::unreachable();
        }
        return result;
    }

    [[using gnu: pure]]
    auto integral(param_type lower_bound, param_type upper_bound) const noexcept -> value_type {
        constexpr param_type kFactor{-(1.0 / 24.0)};
        int sign = 1;
        if (lower_bound > upper_bound) {
            std::swap(lower_bound, upper_bound);
            sign = -1;
        }
        const std::size_t size{m_ts.size()};
        const std::size_t istart =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), lower_bound) - m_ts.cbegin();
        const std::size_t iend =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), upper_bound) - m_ts.cbegin();
        param_type h;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (eval(lower_bound) + eval(upper_bound)) * h * 0.5 * sign;
        }
        value_type result{0.0};
        h = m_ts[istart] - lower_bound;
        result += (eval(lower_bound) + m_ys[istart]) * h * 0.5 +
                  (derivative2(lower_bound) + m_ddys[istart]) * h * h * h * kFactor;
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            h = m_ts[i + 1] - m_ts[i];
            result += (m_ys[i] + m_ys[i + 1]) * h * 0.5 +
                      (m_ddys[i] + m_ddys[i + 1]) * h * h * h * kFactor;
        }
        h = upper_bound - m_ts[iend - 1];
        result += (m_ys[iend - 1] + eval(upper_bound)) * h * 0.5 +
                  (m_ddys[iend - 1] + derivative2(upper_bound)) * h * h * h * kFactor;
        result *= sign;
        return result;
    }

    [[using gnu: pure, always_inline]]
    auto minT() const noexcept -> param_type {
        return m_ts.front();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> param_type {
        return m_ts.back();
    }

    [[using gnu: pure]]
    auto minY() const noexcept -> value_type
        requires std::floating_point<value_type>
    {
        std::size_t pos = std::min_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        param_type h, ratio;
        value_type derivative, derivative2;
        if (pos == 0) {
            h = m_ts[1] - m_ts[0];
            derivative = cuberpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], 0.0, h);
            if (derivative < 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], 1.0, h);
            if (derivative < 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = cuberpd(m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], 0.0, h);
            if (derivative < 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
            derivative2 = lerp(m_ddys[pos - 1], m_ddys[pos], ratio);
            ratio -= derivative / (derivative2 * h);
        }
        return cuberp(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure]]
    auto maxY() const noexcept -> value_type
        requires std::floating_point<value_type>
    {
        std::size_t pos = std::max_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        param_type h, ratio;
        value_type derivative, derivative2;
        if (pos == 0) {
            h = m_ts[1] - m_ts[0];
            derivative = cuberpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], 0.0, h);
            if (derivative > 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], 1.0, h);
            if (derivative > 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = cuberpd(m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], 0.0, h);
            if (derivative > 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
            derivative2 = lerp(m_ddys[pos - 1], m_ddys[pos], ratio);
            ratio -= derivative / (derivative2 * h);
        }
        return cuberp(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(param_type t) const noexcept -> value_type {
        return eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> const std::vector<param_type>& {
        return m_ts;
    }

    [[using gnu: pure, always_inline]]
    auto ys() const noexcept -> const std::vector<value_type>& {
        return m_ys;
    }

    [[using gnu: pure, always_inline]]
    auto ddys() const noexcept -> const std::vector<value_type>& {
        return m_ddys;
    }

  private:
    struct [[nodiscard]] TridiagonalMatrix final {
        [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
        auto luDcmp(std::span<const value_type> b) const noexcept -> std::vector<value_type> {
            const std::size_t mat_size{a_diag.size()};
            std::vector<value_type> x(mat_size);
            std::vector<param_type> u0(mat_size);
            std::vector<param_type> l1(mat_size - 1);
            const std::vector<param_type>& u1{a_up};

            u0[0] = a_diag[0];
            l1[0] = a_low[0] / u0[0];
            for (std::size_t i{1}; i < mat_size - 1; ++i) {
                u0[i] = a_diag[i] - l1[i - 1] * u1[i - 1];
                l1[i] = a_low[i] / u0[i];
            }
            u0[mat_size - 1] = a_diag[mat_size - 1] - l1[mat_size - 2] * u1[mat_size - 2];

            x[0] = b[0];
            for (std::size_t i{1}; i < mat_size; ++i) {
                x[i] = b[i] - l1[i - 1] * x[i - 1];
            }

            x[mat_size - 1] = x[mat_size - 1] / u0[mat_size - 1];
            for (int i = mat_size - 2; i > -1; --i) {
                x[i] = (x[i] - u1[i] * x[i + 1]) / u0[i];
            }
            return x;
        }

        std::vector<param_type> a_low;
        std::vector<param_type> a_diag;
        std::vector<param_type> a_up;
    };

    struct [[nodiscard]] PeriodicTridiagonalMatrix final {
        [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
        auto luDcmp(std::span<const value_type> b) const noexcept -> std::vector<value_type> {
            const std::size_t mat_size{a_diag.size()};
            std::vector<value_type> x(mat_size);
            std::vector<param_type> u0(mat_size);
            std::vector<param_type> l1(mat_size - 1);
            const std::vector<param_type>& u1{a_up};
            std::vector<param_type> l_bottom(mat_size - 2);
            std::vector<param_type> u_top(mat_size - 2);

            u0[0] = a_diag[0];
            l1[0] = a_low[0] / u0[0];
            u_top[0] = a_top;
            l_bottom[0] = a_bottom / u0[0];
            for (std::size_t i{1}; i < mat_size - 2; ++i) {
                u0[i] = a_diag[i] - l1[i - 1] * u1[i - 1];
                l1[i] = a_low[i] / u0[i];
                u_top[i] = -l1[i - 1] * u_top[i - 1];
                l_bottom[i] = -l_bottom[i - 1] * u1[i - 1] / u0[i];
            }
            u0[mat_size - 2] = a_diag[mat_size - 2] - l1[mat_size - 3] * u1[mat_size - 3];
            l1[mat_size - 2] = a_low[mat_size - 2] / u0[mat_size - 2];
            u0[mat_size - 1] = a_diag[mat_size - 1] - l1[mat_size - 2] * u1[mat_size - 2];
            for (std::size_t i{0}; i < mat_size - 2; ++i) {
                u0[mat_size - 1] -= l_bottom[i] * u_top[i];
            }

            x[0] = b[0];
            for (std::size_t i{1}; i < mat_size; ++i) {
                x[i] = b[i] - l1[i - 1] * x[i - 1];
            }
            for (std::size_t i{0}; i < mat_size - 2; ++i) {
                x[mat_size - 1] -= l_bottom[i] * x[i];
            }

            x[mat_size - 1] = x[mat_size - 1] / u0[mat_size - 1];
            x[mat_size - 2] =
                (x[mat_size - 2] - u1[mat_size - 2] * x[mat_size - 1]) / u0[mat_size - 2];
            for (int i = mat_size - 3; i > -1; --i) {
                x[i] = (x[i] - u1[i] * x[i + 1] - u_top[i] * x[mat_size - 1]) / u0[i];
            }
            return x;
        }

        param_type a_bottom;
        std::vector<param_type> a_low;
        std::vector<param_type> a_diag;
        std::vector<param_type> a_up;
        param_type a_top;
    };

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative2(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        const param_type ratio = (t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
        return lerp(m_ddys[pos - 1], m_ddys[pos], ratio);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative3(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        return (m_ddys[pos] - m_ddys[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
    }

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_ts;
        archive & m_ys;
        archive & m_ddys;
        return;
    }

    std::vector<param_type> m_ts{};
    std::vector<value_type> m_ys{};
    std::vector<value_type> m_ddys{};
};

using PiecewiseCubicFunction1f = PiecewiseCubicFunction1<float>;
using PiecewiseCubicFunction1d = PiecewiseCubicFunction1<double>;

} // namespace boyle::math
