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
#include <format>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] PiecewiseLinearFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr U kDuplicateCriterion{kEpsilon};

    [[using gnu: flatten, leaf]]
    explicit PiecewiseLinearFunction1(std::vector<U> ts, std::vector<T> ys) noexcept(
        !BOYLE_CHECK_PARAMS
    ) {
#if BOYLE_CHECK_PARAMS == 1
        if (ts.size() < 2 || ys.size() < 2) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}.",
                ts.size(), ys.size()
            ));
        } else if (ts.size() != ys.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}.",
                ts.size(), ys.size()
            ));
        } else if (!std::ranges::is_sorted(ts.cbegin(), ts.cend())) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        } else if (hasDuplicates(ts.cbegin(), ts.cend(), kDuplicateCriterion)) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }
#endif
        m_ts = std::move(ts);
        m_ys = std::move(ys);
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseLinearFunction1);
    ~PiecewiseLinearFunction1() noexcept = default;

    [[using gnu: pure, always_inline, hot]]
    auto eval(U t) const noexcept -> T {
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            return lerp(m_ys[0], m_ys[1], (t - m_ts[0]) / (m_ts[1] - m_ts[0]));
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.back()
            );
            return lerp(
                m_ys[pos - 1], m_ys[pos - 2], (t - m_ts[pos - 1]) / (m_ts[pos - 2] - m_ts[pos - 1])
            );
        }
        return lerp(m_ys[pos - 1], m_ys[pos], (t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]));
    }

    [[using gnu: pure, always_inline, hot]]
    auto derivative(U t) const noexcept -> T {
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            return (m_ys[1] - m_ys[0]) / (m_ts[1] - m_ts[0]);
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.back()
            );
            return (m_ys[pos - 1] - m_ys[pos - 2]) / (m_ts[pos - 1] - m_ts[pos - 2]);
        }
        return (m_ys[pos] - m_ys[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
    }

    [[using gnu: pure, always_inline]]
    auto derivative(U t, [[maybe_unused]] unsigned int order) const
        noexcept(!BOYLE_CHECK_PARAMS) -> T {
#if BOYLE_CHECK_PARAMS == 1
        if (order != 1) {
            throw std::invalid_argument(std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has first order "
                "derivative: order = {0:d}.",
                order
            ));
        }
#endif
        return derivative(t);
    }

    [[using gnu: pure]]
    auto integral(U lower_bound, U upper_bound) const noexcept -> T {
        int sign{1};
        if (lower_bound > upper_bound) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The lower_bound of integral should always be "
                "less than upper_bound: lower_bound = {0:.6f} while upper_bound = {1:.6f}.",
                lower_bound, upper_bound
            );
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
        T result{0.0};
        result += (eval(lower_bound) + m_ys[istart]) * (m_ts[istart] - lower_bound) * 0.5;
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            result += (m_ys[i] + m_ys[i + 1]) * (m_ts[i + 1] - m_ts[i]) * 0.5;
        }
        result += (m_ys[iend - 1] + eval(upper_bound)) * (upper_bound - m_ts[iend - 1]) * 0.5;
        result *= sign;
        return result;
    }

    [[using gnu: pure, always_inline]]
    auto minT() const noexcept -> U {
        return m_ts.front();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> U {
        return m_ts.back();
    }

    [[using gnu: pure, always_inline]]
    auto minY() const noexcept -> T
        requires std::floating_point<T>
    {
        return *std::min_element(m_ys.cbegin(), m_ys.cend());
    }

    [[using gnu: pure, always_inline]]
    auto maxY() const noexcept -> T
        requires std::floating_point<T>
    {
        return *std::max_element(m_ys.cbegin(), m_ys.cend());
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto eval(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<T>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<T> ys;
        ys.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            ys.emplace_back(eval(*it));
        }
        return ys;
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto derivative(ForwardIt first, std::sentinel_for<ForwardIt> auto last) const noexcept
        -> std::vector<T>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<T> dys;
        dys.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            dys.emplace_back(derivative(*it));
        }
        return dys;
    }

    template <std::forward_iterator ForwardIt>
    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto derivative(ForwardIt first, std::sentinel_for<ForwardIt> auto last, unsigned int order)
        const -> std::vector<T>
        requires std::same_as<T, typename std::iterator_traits<ForwardIt>::value_type>
    {
        std::vector<T> dnys;
        dnys.reserve(last - first);
        for (ForwardIt it{first}; it != last; ++it) {
            dnys.emplace_back(derivative(*it, order));
        }
        return dnys;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(U t) const noexcept -> T {
        return eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> const std::vector<U>& {
        return m_ts;
    }

    [[using gnu: pure, always_inline]]
    auto ys() const noexcept -> const std::vector<T>& {
        return m_ys;
    }

  private:
    std::vector<U> m_ts{};
    std::vector<T> m_ys{};

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
