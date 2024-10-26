/**
 * @file piecewise_linear_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-10
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <ranges>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "fmt/format.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"
#include "boyle/math/vec3.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] PiecewiseLinearFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr param_type kDuplicateCriterion{1E-8};

    PiecewiseLinearFunction1() noexcept = default;
    PiecewiseLinearFunction1(const PiecewiseLinearFunction1& other) noexcept = default;
    auto operator=(const PiecewiseLinearFunction1& other
    ) noexcept -> PiecewiseLinearFunction1& = default;
    PiecewiseLinearFunction1(PiecewiseLinearFunction1&& other) noexcept = default;
    auto operator=(PiecewiseLinearFunction1&& other
    ) noexcept -> PiecewiseLinearFunction1& = default;
    ~PiecewiseLinearFunction1() noexcept = default;

    [[using gnu: ]]
    explicit PiecewiseLinearFunction1(
        std::vector<param_type> ts, std::vector<value_type> ys
    ) noexcept(!BOYLE_CHECK_PARAMS)
        : m_ts{std::move(ts)}, m_ys{std::move(ys)} {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ts.size() < 2 || m_ys.size() < 2) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}.",
                m_ts.size(), m_ys.size()
            ));
        }
        if (m_ts.size() != m_ys.size()) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}.",
                m_ts.size(), m_ys.size()
            ));
        }
        if (!std::ranges::is_sorted(m_ts.cbegin(), m_ts.cend())) [[unlikely]] {
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
#endif
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto eval(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return lerp(m_ys[0], m_ys[1], (t - m_ts[0]) / (m_ts[1] - m_ts[0]));
        }
        if (pos == m_ts.size()) {
            return lerp(
                m_ys[pos - 1], m_ys[pos - 2], (t - m_ts[pos - 1]) / (m_ts[pos - 2] - m_ts[pos - 1])
            );
        }
        return lerp(m_ys[pos - 1], m_ys[pos], (t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]));
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return (m_ys[1] - m_ys[0]) / (m_ts[1] - m_ts[0]);
        }
        if (pos == m_ts.size()) {
            return (m_ys[pos - 1] - m_ys[pos - 2]) / (m_ts[pos - 1] - m_ts[pos - 2]);
        }
        return (m_ys[pos] - m_ys[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
    }

    [[using gnu: pure, always_inline]]
    auto derivative(param_type t, [[maybe_unused]] unsigned int order) const
        noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (order != 1) [[unlikely]] {
            throw std::invalid_argument(fmt::format(
                "Invalid argument error! The PiecewiseLinearFunction only has first order "
                "derivative: order = {0:d}.",
                order
            ));
        }
#endif
        return derivative(t);
    }

    [[using gnu: pure]]
    auto integral(param_type lower_bound, param_type upper_bound) const noexcept -> value_type {
        int sign{1};
        if (lower_bound > upper_bound) {
            std::swap(lower_bound, upper_bound);
            sign = -1;
        }
        const std::size_t size{m_ts.size()};
        const std::size_t istart =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), lower_bound) - m_ts.cbegin();
        const std::size_t iend =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), upper_bound) - m_ts.cbegin();
        if (istart == size || iend == 0 || istart == iend) {
            return (eval(lower_bound) + eval(upper_bound)) * (upper_bound - lower_bound) * 0.5 *
                   sign;
        }
        value_type result{0.0};
        result += (eval(lower_bound) + m_ys[istart]) * (m_ts[istart] - lower_bound) * 0.5;
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            result += (m_ys[i] + m_ys[i + 1]) * (m_ts[i + 1] - m_ts[i]) * 0.5;
        }
        result += (m_ys[iend - 1] + eval(upper_bound)) * (upper_bound - m_ts[iend - 1]) * 0.5;
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

    [[using gnu: pure, always_inline]]
    auto minY() const noexcept -> value_type
        requires std::floating_point<value_type>
    {
        return *std::min_element(m_ys.cbegin(), m_ys.cend());
    }

    [[using gnu: pure, always_inline]]
    auto maxY() const noexcept -> value_type
        requires std::floating_point<value_type>
    {
        return *std::max_element(m_ys.cbegin(), m_ys.cend());
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

  private:
    std::vector<param_type> m_ts{};
    std::vector<value_type> m_ys{};

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_ts;
        archive & m_ys;
        return;
    }
};

using PiecewiseLinearFunction1f = PiecewiseLinearFunction1<float>;
using PiecewiseLinearFunction1d = PiecewiseLinearFunction1<double>;

} // namespace boyle::math
