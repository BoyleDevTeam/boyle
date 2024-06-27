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
#include <format>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/logging.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/concepts.hpp"
#include "boyle/math/cubic_interpolation.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] PiecewiseCubicFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr U kDuplicateCriterion{kEpsilon};

    struct BoundaryMode final {
        unsigned int order;
        T derivative;
    };

    [[using gnu: always_inline]]
    explicit PiecewiseCubicFunction1(std::vector<U> ts, std::vector<T> ys)
        : PiecewiseCubicFunction1<T, U>(
              std::move(ts), std::move(ys), BoundaryMode{2, 0.0}, BoundaryMode{2, 0.0}
          ) {}

    [[using gnu: flatten, leaf]]
    explicit PiecewiseCubicFunction1(
        std::vector<U> ts, std::vector<T> ys, BoundaryMode b0, BoundaryMode bf
    ) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
        if (ts.size() < 2 || ys.size() < 2) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                ts.size(), ys.size()
            ));
        } else if (ts.size() != ys.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                ts.size(), ys.size()
            ));
        } else if (!std::is_sorted(ts.cbegin(), ts.cend())) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        } else if (hasDuplicates(ts.cbegin(), ts.cend(), kDuplicateCriterion)) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }

        if (b0.order == 0 || b0.order > 2) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "b0.order = {0:d}.",
                b0.order
            ));
        }
        if (bf.order == 0 || bf.order > 2) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "bf.order = {0:d}.",
                bf.order
            ));
        }
#endif

        const std::size_t size{ts.size()};
        std::vector<U> hs(size - 1);
        std::vector<T> ds(size - 1);
        std::vector<U> a_diag(size);
        std::vector<T> b(size);

        hs[0] = ts[1] - ts[0];
        ds[0] = (ys[1] - ys[0]) / hs[0];
        for (std::size_t i{1}; i < size - 1; i++) {
            hs[i] = ts[i + 1] - ts[i];
            ds[i] = (ys[i + 1] - ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }

        std::vector<U> a_low{hs};
        std::vector<U> a_up{hs};

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

        m_ts = std::move(ts);
        m_ys = std::move(ys);
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseCubicFunction1);
    ~PiecewiseCubicFunction1() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto eval(U t) const noexcept -> T {
        constexpr std::array<U, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            const U h{m_ts[1] - m_ts[0]};
            const U ratio{(t - m_ts[0]) / h};
            return lerp(m_ys[0], m_ys[1], ratio) +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * (t - m_ts[0]) * h;
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, m_ts.back()
            );
            const U h{m_ts[pos - 2] - m_ts[pos - 1]};
            const U ratio{(t - m_ts[pos - 1]) / h};
            return lerp(m_ys[pos - 1], m_ys[pos - 2], ratio) +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) *
                       (t - m_ts[pos - 1]) * h;
        }
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return cuberp(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    auto derivative(U t) const noexcept -> T {
        constexpr std::array<U, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            const U h{m_ts[1] - m_ts[0]};
            return (m_ys[1] - m_ys[0]) / h +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * h;
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, m_ts.back()
            );
            const U h{m_ts[pos - 2] - m_ts[pos - 1]};
            return (m_ys[pos - 2] - m_ys[pos - 1]) / h +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) * h;
        }
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio{(t - m_ts[pos - 1]) / h};
        return cuberpd(m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    auto derivative(U t, unsigned int order) const noexcept(!BOYLE_CHECK_PARAMS) -> T {
#if BOYLE_CHECK_PARAMS == 1
        if (order < 1 || order > 3) {
            throw std::invalid_argument(std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has 1, 2, 3 order "
                "derivatives: order = {0:d}.",
                order
            ));
        }
#endif
        T result;
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
        }
        return result;
    }

    [[using gnu: pure]]
    auto integral(U lower_bound, U upper_bound) const noexcept -> T {
        constexpr U kFactor{-(1.0 / 24.0)};
        int sign = 1;
        if (lower_bound > upper_bound) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! the lower_bound of integral should always be "
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
        U h;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (eval(lower_bound) + eval(upper_bound)) * h * 0.5 * sign;
        }
        T result{0.0};
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
    auto minT() const noexcept -> U {
        return m_ts.front();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> U {
        return m_ts.back();
    }

    [[using gnu: pure, flatten, leaf]]
    auto minY() const noexcept -> T
        requires std::floating_point<T>
    {
        std::size_t pos = std::min_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        U h, ratio;
        T derivative, derivative2;
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

    [[using gnu: pure, flatten, leaf]]
    auto maxY() const noexcept -> T
        requires std::floating_point<T>
    {
        std::size_t pos = std::max_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        U h, ratio;
        T derivative, derivative2;
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

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> const std::vector<U>& {
        return m_ts;
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(U t) const noexcept -> T {
        return eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto ys() const noexcept -> const std::vector<T>& {
        return m_ys;
    }

    [[using gnu: pure, always_inline]]
    auto ddys() const noexcept -> const std::vector<T>& {
        return m_ddys;
    }

  private:
    struct [[nodiscard]] TridiagonalMatrix final {
        [[using gnu: pure, always_inline]] [[nodiscard]]
        auto luDcmp(const std::vector<T>& b) const noexcept -> std::vector<T> {
            const std::size_t mat_size{a_diag.size()};
            std::vector<T> x(mat_size);
            std::vector<U> u0(mat_size);
            std::vector<U> l1(mat_size - 1);
            const std::vector<U>& u1{a_up};

            u0[0] = a_diag[0];
            l1[0] = a_low[0] / u0[0];
            for (std::size_t i{1}; i < mat_size - 1; i++) {
                u0[i] = a_diag[i] - l1[i - 1] * u1[i - 1];
                l1[i] = a_low[i] / u0[i];
            }
            u0[mat_size - 1] = a_diag[mat_size - 1] - l1[mat_size - 2] * u1[mat_size - 2];

            x[0] = b[0];
            for (std::size_t i{1}; i < mat_size; i++) {
                x[i] = b[i] - l1[i - 1] * x[i - 1];
            }

            x[mat_size - 1] = x[mat_size - 1] / u0[mat_size - 1];
            for (std::int64_t i = mat_size - 2; i > -1; i--) {
                x[i] = (x[i] - u1[i] * x[i + 1]) / u0[i];
            }
            return x;
        }

        std::vector<U> a_low;
        std::vector<U> a_diag;
        std::vector<U> a_up;
    };

    [[using gnu: flatten, leaf]]
    explicit PiecewiseCubicFunction1(
        std::vector<U> ts, std::vector<T> ys, std::vector<T> ddys
    ) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
        if (ts.size() < 2 || ys.size() < 2) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                ts.size(), ys.size()
            ));
        } else if (ts.size() != ys.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                ts.size(), ys.size()
            ));
        } else if (!std::is_sorted(ts.cbegin(), ts.cend())) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        } else if (hasDuplicates(ts.cbegin(), ts.cend(), kDuplicateCriterion)) {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }
        if (ddys.size() != ts.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! ts, ys, ddys must share the same size: ts.size() = "
                "{0:d} while ddys.size() = {1:d}",
                ts.size(), ddys.size()
            ));
        }
#endif
        m_ts = std::move(ts);
        m_ys = std::move(ys);
        m_ddys = std::move(ddys);
    }

    [[using gnu: pure, always_inline]]
    auto derivative2(U t) const noexcept -> T {
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            return T{0.0};
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, m_ts.back()
            );
            return T{0.0};
        }
        const U ratio = (t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
        return lerp(m_ddys[pos - 1], m_ddys[pos], ratio);
    }

    [[using gnu: pure, always_inline]]
    auto derivative3(U t) const noexcept -> T {
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            return T{0.0};
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, m_ts.back()
            );
            return T{0.0};
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

    std::vector<U> m_ts{};
    std::vector<T> m_ys{};
    std::vector<T> m_ddys{};
};

using PiecewiseCubicFunction1f = PiecewiseCubicFunction1<float>;
using PiecewiseCubicFunction1d = PiecewiseCubicFunction1<double>;

} // namespace boyle::math
