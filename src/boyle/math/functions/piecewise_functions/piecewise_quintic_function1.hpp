/**
 * @file piecewise_quintic_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-18
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
#include "boyle/math/quintic_interpolation.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class PiecewiseQuinticFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr U kDuplicateCriterion{kEpsilon};

    struct [[nodiscard]] BoundaryMode final {
        unsigned int order;
        T derivative;
    };

    [[using gnu: always_inline]]
    explicit PiecewiseQuinticFunction1(std::vector<U> ts, std::vector<T> ys)
        : PiecewiseQuinticFunction1(
              std::move(ts), std::move(ys),
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}},
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}}
          ) {}

    [[using gnu: flatten, leaf]]
    explicit PiecewiseQuinticFunction1(
        std::vector<U> ts, std::vector<T> ys, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf
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

        if (b0[0].order == 0 || b0[0].order > 4) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: b01.order = {0:d}.",
                b0[0].order
            ));
        }
        if (b0[1].order == 0 || b0[1].order > 4) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: b02.order = {0:d}.",
                b0[0].order
            ));
        }
        if (b0[0].order == b0[1].order) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! You should not set the same order of boundary mode "
                "twice: b01.order == {0:d} while b02.order == {1:d}",
                b0[0].order, b0[1].order
            ));
        }
        if (bf[0].order == 0 || bf[0].order > 4) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: bf1.order = {0:d}.",
                bf[0].order
            ));
        }
        if (bf[1].order == 0 || bf[1].order > 4) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: bf2.order = {0:d}.",
                bf[1].order
            ));
        }
        if (bf[0].order == bf[1].order) {
            throw std::invalid_argument(std::format(
                "Invalid argument detected! You should not set the same order of boundary mode "
                "twice: bf1.order == {0:d} while bf2.order == {1:d}",
                bf[0].order, bf[1].order
            ));
        }
#endif
        constexpr std::array<U, 3> kFactors{-(7.0 / 60.0), -(2.0 / 15.0), -(7.0 / 60.0)};

        if (b0[0].order > b0[1].order) {
            std::swap(b0[0], b0[1]);
        }
        if (bf[0].order > bf[1].order) {
            std::swap(bf[0], bf[1]);
        }

        const std::size_t size{ts.size()};
        std::vector<U> hs(size - 1);
        std::vector<T> ds(size - 1);
        std::vector<U> a_diag(size * 2);
        std::vector<U> a_low_1(size * 2 - 1);
        std::vector<U> a_low_2(size + 1);
        std::vector<U> a_low_3(size);
        std::vector<U> a_low_4(size - 1);
        std::vector<U> a_up_1(size * 2 - 1);
        std::vector<U> a_up_2(size + 1);
        std::vector<U> a_up_3(size);
        std::vector<U> a_up_4(size - 1);
        std::vector<T> b(size * 2);

        hs[0] = ts[1] - ts[0];
        ds[0] = (ys[1] - ys[0]) / hs[0];
        a_low_2[0] = 0.0;
        a_up_2[0] = 0.0;
        for (std::size_t i{1}; i < size - 1; ++i) {
            hs[i] = ts[i + 1] - ts[i];
            ds[i] = (ys[i + 1] - ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            a_low_2[i + 1] = -(6.0 / hs[i]);
            a_low_3[i] = 6.0 / hs[i] + 6.0 / hs[i - 1];
            a_low_4[i - 1] = -(6.0 / hs[i - 1]);
            a_up_2[i] = hs[i - 1] * hs[i - 1] * hs[i - 1] * kFactors[0];
            a_up_3[i] = (hs[i - 1] * hs[i - 1] * hs[i - 1] + hs[i] * hs[i] * hs[i]) * kFactors[1];
            a_up_4[i] = hs[i] * hs[i] * hs[i] * kFactors[2];
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }
        a_low_2[size] = 0.0;
        a_up_2[size] = 0.0;

        std::copy(a_diag.cbegin(), a_diag.cbegin() + size, a_diag.begin() + size);
        std::copy(hs.cbegin(), hs.cend(), a_low_1.begin());
        std::copy(hs.cbegin(), hs.cend(), a_low_1.begin() + size);
        std::copy(hs.cbegin(), hs.cend(), a_up_1.begin());
        std::copy(hs.cbegin(), hs.cend(), a_up_1.begin() + size);
        std::fill(b.begin() + size, b.end(), T{0.0});

        a_low_1[size - 1] = 0.0;
        a_up_1[size - 1] = 0.0;

        if (b0[0].order == 1 && b0[1].order == 2) {
            const U h2 = hs[0] * hs[0];
            a_diag[size] = 1.0;
            a_low_2[1] = -7.5 / h2;
            a_low_3[0] = -15.0 / h2;
            a_up_1[size] = 0.875;
            b[size] = (b0[0].derivative - ds[0]) / (h2 * hs[0]) * 45.0;
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[1].derivative;
        } else if (b0[0].order == 2 && b0[1].order == 4) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[0].derivative;
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        } else if (b0[0].order == 1 && b0[1].order == 3) {
            const U h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[0].derivative) / hs[0] * 3.0;
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / h2);
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[1].derivative / hs[0] * 3.0);
        } else if (b0[0].order == 2 && b0[1].order == 3) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[0].derivative;
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / (hs[0] * hs[0]));
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[1].derivative / hs[0] * 3.0);
        } else if (b0[0].order == 3 && b0[1].order == 4) {
            const U h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = -1.0;
            a_up_3[0] = h2 / 3.0;
            a_up_4[0] = h2 / 6.0;
            b[0] = -b0[0].derivative * hs[0];
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        } else {
            const U h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[0].derivative) / hs[0] * 3.0;
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        }

        if (bf[0].order == 1 && bf[1].order == 2) {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.875;
            a_low_3[size - 1] = -15.0 / h2;
            a_low_4[size - 2] = -7.5 / h2;
            b[size * 2 - 1] = (ds[size - 2] - bf[0].derivative) / (h2 * hs[size - 2]) * 45.0;
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[1].derivative;
        } else if (bf[0].order == 2 && bf[1].order == 4) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[0].derivative;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        } else if (bf[0].order == 1 && bf[1].order == 3) {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[0].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / h2;
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[1].derivative / hs[size - 2] * 3.0;
        } else if (bf[0].order == 2 && bf[1].order == 3) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[0].derivative;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / (hs[size - 2] * hs[size - 2]);
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[1].derivative / hs[size - 2] * 3.0;
        } else if (bf[0].order == 3 && bf[1].order == 4) {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = -1.0;
            a_up_2[size - 1] = h2 / 6.0;
            a_up_3[size - 1] = h2 / 3.0;
            b[size - 1] = bf[0].derivative * hs[size - 2];
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        } else {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[0].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        }

        const OutriggerMatrix A{std::move(a_low_4), std::move(a_low_3), std::move(a_low_2),
                                std::move(a_low_1), std::move(a_diag),  std::move(a_up_1),
                                std::move(a_up_2),  std::move(a_up_3),  std::move(a_up_4)};

        std::vector<T> x = A.gaussSeidel(b);

        m_ddys.resize(size);
        m_d4ys.resize(size);
        std::copy(x.cbegin(), x.cbegin() + size, m_ddys.begin());
        std::copy(x.cbegin() + size, x.cend(), m_d4ys.begin());

        m_ts = std::move(ts);
        m_ys = std::move(ys);
    }

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseQuinticFunction1);
    ~PiecewiseQuinticFunction1() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto eval(U t) const noexcept -> T {
        constexpr std::array<U, 4> kFactors{-(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0};
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
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * (t - m_ts[0]) * h +
                   (m_d4ys[0] * kFactors[2] + m_d4ys[1] * kFactors[3]) * (t - m_ts[0]) * h * h * h;
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
                       (t - m_ts[pos - 1]) * h +
                   (m_d4ys[pos - 1] * kFactors[2] + m_d4ys[pos - 2] * kFactors[3]) *
                       (t - m_ts[pos - 1]) * h * h * h;
        }
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    auto derivative(U t) const noexcept -> T {
        constexpr std::array<U, 4> kFactors{-(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0};
        const std::size_t pos = nearestUpperElement(m_ts.cbegin(), m_ts.cend(), t) - m_ts.cbegin();
        if (pos == 0) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, m_ts.front()
            );
            const U h{m_ts[1] - m_ts[0]};
            return (m_ys[1] - m_ys[0]) / h +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * h +
                   (m_d4ys[0] * kFactors[2] + m_d4ys[1] * kFactors[3]) * h * h * h;
        }
        if (pos == m_ts.size()) {
            BOYLE_LOG_WARN(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, m_ts.back()
            );
            const U h{m_ts[pos - 2] - m_ts[pos - 1]};
            return (m_ys[pos - 2] - m_ys[pos - 1]) / h +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) * h +
                   (m_d4ys[pos - 1] * kFactors[2] + m_d4ys[pos - 2] * kFactors[3]) * h * h * h;
        }
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio{(t - m_ts[pos - 1]) / h};
        return quinerpd(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    auto derivative(U t, unsigned int order) const noexcept(!BOYLE_CHECK_PARAMS) -> T {
#if BOYLE_CHECK_PARAMS == 1
        if (order < 1 || order > 5) {
            throw std::invalid_argument(std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has 1st, 2nd, 3rd, 4th "
                "order derivatives: order = {0:d}.",
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
        case 4:
            result = derivative4(t);
            break;
        case 5:
            result = derivative5(t);
            break;
        }
        return result;
    }

    [[using gnu: pure, flatten, leaf]]
    auto integral(U lower_bound, U upper_bound) const noexcept -> T {
        constexpr std::array<U, 2> kFactors{-(1.0 / 24.0), 1.0 / 240.0};
        int sign{1};
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
        U h, h3;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (eval(lower_bound) + eval(upper_bound)) * h * 0.5 * sign;
        }
        T result{0.0};
        h = m_ts[istart] - lower_bound;
        h3 = h * h * h;
        result += (eval(lower_bound) + m_ys[istart]) * h * 0.5 +
                  (derivative2(lower_bound) + m_ddys[istart]) * h3 * kFactors[0] +
                  (derivative4(lower_bound) + m_d4ys[istart]) * h3 * h * h * kFactors[1];
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            h = m_ts[i + 1] - m_ts[i];
            h3 = h * h * h;
            result += (m_ys[i] + m_ys[i + 1]) * h * 0.5 +
                      (m_ddys[i] + m_ddys[i + 1]) * h3 * kFactors[0] +
                      (m_d4ys[i] + m_d4ys[i + 1]) * h3 * h * h * kFactors[1];
        }
        h = upper_bound - m_ts[iend - 1];
        h3 = h * h * h;
        result += (m_ys[iend - 1] + eval(upper_bound)) * h * 0.5 +
                  (m_ddys[iend - 1] + derivative2(upper_bound)) * h3 * kFactors[0] +
                  (m_d4ys[iend - 1] + derivative4(upper_bound)) * h3 * h * h * kFactors[1];
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
            derivative = quinerpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], m_ys[0], m_ys[1], 0.0, h);
            if (derivative < 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], 1.0, h
            );
            if (derivative < 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = quinerpd(
                m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], m_d4ys[pos],
                m_d4ys[pos + 1], 0.0, h
            );
            if (derivative < 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], ratio, h
            );
            derivative2 =
                cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
            ratio = derivative / (derivative2 * h);
        }
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
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
            derivative = quinerpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], m_ys[0], m_ys[1], 0.0, h);
            if (derivative > 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], 1.0, h
            );
            if (derivative > 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = quinerpd(
                m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], m_d4ys[pos],
                m_d4ys[pos + 1], 0.0, h
            );
            if (derivative > 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], ratio, h
            );
            derivative2 =
                cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
            ratio = derivative / (derivative2 * h);
        }
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
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

    [[using gnu: pure, always_inline]]
    auto ddys() const noexcept -> const std::vector<T>& {
        return m_ddys;
    }

    [[using gnu: pure, always_inline]]
    auto d4ys() const noexcept -> const std::vector<T>& {
        return m_d4ys;
    }

  private:
    struct [[nodiscard]] OutriggerMatrix final {
        [[using gnu: pure, always_inline]] [[nodiscard]]
        auto gaussSeidel(const std::vector<T>& b) const noexcept -> std::vector<T> {
            const std::size_t mat_size{a_diag.size()};
            const std::size_t half_mat_size{mat_size / 2};
            std::vector<T> x(mat_size, 0.0);
            for (std::size_t num_iter{40}; num_iter != 0U; --num_iter) {
                x[0] = (b[0] - a_up_1[0] * x[1] - a_up_2[0] * x[half_mat_size - 1] -
                        a_up_3[0] * x[half_mat_size] - a_up_4[0] * x[half_mat_size + 1]) /
                       a_diag[0];
                for (std::size_t i{1}; i < half_mat_size - 1; ++i) {
                    x[i] =
                        (b[i] - a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1] -
                         a_up_2[i] * x[half_mat_size + i - 1] - a_up_3[i] * x[half_mat_size + i] -
                         a_up_4[i] * x[half_mat_size + i + 1]) /
                        a_diag[i];
                }
                x[half_mat_size - 1] = (b[half_mat_size - 1] - a_low_2[0] * x[0] -
                                        a_low_1[half_mat_size - 2] * x[half_mat_size - 2] -
                                        a_up_1[half_mat_size - 1] * x[half_mat_size] -
                                        a_up_2[half_mat_size - 1] * x[mat_size - 2] -
                                        a_up_3[half_mat_size - 1] * x[mat_size - 1]) /
                                       a_diag[half_mat_size - 1];
                x[half_mat_size] = (b[half_mat_size] - a_low_3[0] * x[0] - a_low_2[1] * x[1] -
                                    a_low_1[half_mat_size - 1] * x[half_mat_size - 1] -
                                    a_up_1[half_mat_size] * x[half_mat_size + 1] -
                                    a_up_2[half_mat_size] * x[mat_size - 1]) /
                                   a_diag[half_mat_size];
                for (std::size_t i{half_mat_size + 1}; i < mat_size - 1; ++i) {
                    x[i] = (b[i] - a_low_4[i - half_mat_size - 1] * x[i - half_mat_size - 1] -
                            a_low_3[i - half_mat_size] * x[i - half_mat_size] -
                            a_low_2[i - half_mat_size + 1] * x[i - half_mat_size + 1] -
                            a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1]) /
                           a_diag[i];
                }
                x[mat_size - 1] =
                    (b[mat_size - 1] - a_low_4[half_mat_size - 2] * x[half_mat_size - 2] -
                     a_low_3[half_mat_size - 1] * x[half_mat_size - 1] -
                     a_low_2[half_mat_size] * x[half_mat_size] -
                     a_low_1[mat_size - 2] * x[mat_size - 2]) /
                    a_diag[mat_size - 1];
            }
            return x;
        }

        std::vector<U> a_low_4;
        std::vector<U> a_low_3;
        std::vector<U> a_low_2;
        std::vector<U> a_low_1;
        std::vector<U> a_diag;
        std::vector<U> a_up_1;
        std::vector<U> a_up_2;
        std::vector<U> a_up_3;
        std::vector<U> a_up_4;
    };

    [[using gnu: flatten, leaf]]
    explicit PiecewiseQuinticFunction1(
        std::vector<U> ts, std::vector<T> ys, std::vector<T> ddys, std::vector<T> d4ys
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
        if (d4ys.size() != ts.size()) {
            throw std::invalid_argument(std::format(
                "Invalid arguments detected! ts, ys, ddys must share the same size: ts.size() = "
                "{0:d} while d4ys.size() = {1:d}",
                ts.size(), d4ys.size()
            ));
        }
#endif
        m_ts = std::move(ts);
        m_ys = std::move(ys);
        m_ddys = std::move(ddys);
        m_d4ys = std::move(d4ys);
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
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio = (t - m_ts[pos - 1]) / h;
        return cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
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
        const U h{m_ts[pos] - m_ts[pos - 1]};
        const U ratio{(t - m_ts[pos - 1]) / h};
        return cuberpd(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    auto derivative4(U t) const noexcept -> T {
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
        const U ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return lerp(m_d4ys[pos - 1], m_d4ys[pos], ratio);
    }

    [[using gnu: pure, always_inline]]
    auto derivative5(U t) const noexcept -> T {
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
        return (m_d4ys[pos] - m_d4ys[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
    }

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_ts;
        archive & m_ys;
        archive & m_ddys;
        archive & m_d4ys;
        return;
    }

    std::vector<U> m_ts{};
    std::vector<T> m_ys{};
    std::vector<T> m_ddys{};
    std::vector<T> m_d4ys{};
};

using PiecewiseQuinticFunction1f = PiecewiseQuinticFunction1<float>;
using PiecewiseQuinticFunction1d = PiecewiseQuinticFunction1<double>;

} // namespace boyle::math
