/**
 * @file piecewise_quintic_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <format>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "common/utils/macros.hpp"
#include "math/cubic_interpolation.hpp"
#include "math/functions/function1.hpp"
#include "math/quintic_interpolation.hpp"
#include "math/type_traits.hpp"
#include "math/utils.hpp"

namespace tiny_pnc {
namespace math {

template <typename T, typename U = T>
class PiecewiseQuinticFunction1 final : public Function1<PiecewiseQuinticFunction1<T, U>, T, U> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    static constexpr U kDuplicateCriterion{kEpsilon};

    struct [[nodiscard]] BoundaryMode final {
        unsigned int order;
        T derivative;
    };

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseQuinticFunction1);

    [[using gnu: always_inline]] explicit PiecewiseQuinticFunction1(
        std::vector<U> ts, std::vector<T> ys
    )
        : PiecewiseQuinticFunction1(
              std::move(ts), std::move(ys),
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}},
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}}
          ) {}

    [[using gnu: flatten, leaf]] explicit PiecewiseQuinticFunction1(
        std::vector<U> ts, std::vector<T> ys, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf
    ) {
        constexpr std::array<U, 3> kFactors{-(7.0 / 60.0), -(2.0 / 15.0), -(7.0 / 60.0)};
        if (ts.size() < 2 || ys.size() < 2) {
            std::string error_msg = std::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                ts.size(), ys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        } else if (ts.size() != ys.size()) {
            std::string error_msg = std::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                ts.size(), ys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        } else if (!std::is_sorted(ts.cbegin(), ts.cend())) {
            std::string error_msg =
                std::format("Invalid arguments detected! ts has to be a sorted array!");
            throw std::invalid_argument(std::move(error_msg));
        } else if (hasDuplicates(ts.cbegin(), ts.cend(), kDuplicateCriterion)) {
            std::string error_msg =
                std::format("Invalid arguments detected! ts can not have duplicated elements!");
            throw std::invalid_argument(std::move(error_msg));
        }

        ts_ = std::move(ts);
        ys_ = std::move(ys);

        if (b0[0].order == 0 || b0[0].order > 4) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: b01.order = {0:d}.",
                b0[0].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (b0[1].order == 0 || b0[1].order > 4) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: b02.order = {0:d}.",
                b0[0].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (b0[0].order == b0[1].order) {
            std::string error_msg = std::format(
                "Invalid argument detected! You should not set the same order of boundary mode "
                "twice: b01.order == {0:d} while b02.order == {1:d}",
                b0[0].order, b0[1].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (bf[0].order == 0 || bf[0].order > 4) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: bf1.order = {0:d}.",
                bf[0].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (bf[1].order == 0 || bf[1].order > 4) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be either 1, 2, 3, "
                "4, 5: bf2.order = {0:d}.",
                bf[1].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (bf[0].order == bf[1].order) {
            std::string error_msg = std::format(
                "Invalid argument detected! You should not set the same order of boundary mode "
                "twice: bf1.order == {0:d} while bf2.order == {1:d}",
                bf[0].order, bf[1].order
            );
            throw std::invalid_argument(std::move(error_msg));
        }

        const std::size_t size{ts_.size()};
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

        hs[0] = ts_[1] - ts_[0];
        ds[0] = (ys_[1] - ys_[0]) / hs[0];
        a_low_2[0] = 0.0;
        a_up_2[0] = 0.0;
        for (std::size_t i{1}; i < size - 1; ++i) {
            hs[i] = ts_[i + 1] - ts_[i];
            ds[i] = (ys_[i + 1] - ys_[i]) / hs[i];
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

        if (b0[0].order == 2) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[0].derivative;
        } else if (b0[0].order == 4) {
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[0].derivative;
        } else if (b0[0].order == 1) {
            const U h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[0].derivative) / hs[0] * 3.0;
        } else {
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / (hs[0] * hs[0]));
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[0].derivative / hs[0] * 3.0);
        }

        if (b0[1].order == 2) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[1].derivative;
        } else if (b0[1].order == 4) {
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        } else if (b0[1].order == 1) {
            const U h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[1].derivative) / hs[0] * 3.0;
        } else {
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / (hs[0] * hs[0]));
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[1].derivative / hs[0] * 3.0);
        }

        if (bf[0].order == 2) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[0].derivative;
        } else if (bf[0].order == 4) {
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[0].derivative;
        } else if (bf[0].order == 1) {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[0].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
        } else {
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / (hs[size - 2] * hs[size - 2]);
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[0].derivative / hs[size - 2] * 3.0;
        }

        if (bf[1].order == 2) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[1].derivative;
        } else if (bf[1].order == 4) {
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        } else if (bf[1].order == 1) {
            const U h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[1].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
        } else {
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / (hs[size - 2] * hs[size - 2]);
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[1].derivative / hs[size - 2] * 3.0;
        }

        const OutriggerMatrix A{std::move(a_low_4), std::move(a_low_3), std::move(a_low_2),
                                std::move(a_low_1), std::move(a_diag),  std::move(a_up_1),
                                std::move(a_up_2),  std::move(a_up_3),  std::move(a_up_4)};

        std::vector<T> x = A.gaussSeidel(b);

        ddys_.resize(size);
        d4ys_.resize(size);
        std::copy(x.cbegin(), x.cbegin() + size, ddys_.begin());
        std::copy(x.cbegin() + size, x.cend(), d4ys_.begin());
    }

    ~PiecewiseQuinticFunction1() noexcept override = default;

    [[using gnu: pure, always_inline]] T operator()(U t) const noexcept {
        constexpr std::array<U, 4> kFactors{-(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0};
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            const U h{ts_[1] - ts_[0]};
            const U ratio{(t - ts_[0]) / h};
            return lerp(ys_[0], ys_[1], ratio) +
                   (ddys_[0] * kFactors[0] + ddys_[1] * kFactors[1]) * (t - ts_[0]) * h +
                   (d4ys_[0] * kFactors[2] + d4ys_[1] * kFactors[3]) * (t - ts_[0]) * h * h * h;
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            const U h{ts_[pos - 2] - ts_[pos - 1]};
            const U ratio{(t - ts_[pos - 1]) / h};
            return lerp(ys_[pos - 1], ys_[pos - 2], ratio) +
                   (ddys_[pos - 1] * kFactors[0] + ddys_[pos - 2] * kFactors[1]) *
                       (t - ts_[pos - 1]) * h +
                   (d4ys_[pos - 1] * kFactors[2] + d4ys_[pos - 2] * kFactors[3]) *
                       (t - ts_[pos - 1]) * h * h * h;
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio{(t - ts_[pos - 1]) / (ts_[pos] - ts_[pos - 1])};
        return quinerp(
            ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t) const noexcept {
        constexpr std::array<U, 4> kFactors{-(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0};
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            const U h{ts_[1] - ts_[0]};
            return (ys_[1] - ys_[0]) / h + (ddys_[0] * kFactors[0] + ddys_[1] * kFactors[1]) * h +
                   (d4ys_[0] * kFactors[2] + d4ys_[1] * kFactors[3]) * h * h * h;
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            const U h{ts_[pos - 2] - ts_[pos - 1]};
            return (ys_[pos - 2] - ys_[pos - 1]) / h +
                   (ddys_[pos - 1] * kFactors[0] + ddys_[pos - 2] * kFactors[1]) * h +
                   (d4ys_[pos - 1] * kFactors[2] + d4ys_[pos - 2] * kFactors[3]) * h * h * h;
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio{(t - ts_[pos - 1]) / h};
        return quinerpd(
            ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t, unsigned int order) const {
        T result;
        if (order == 1) {
            result = derivative(t);
        } else if (order == 2) {
            result = derivative2(t);
        } else if (order == 3) {
            result = derivative3(t);
        } else if (order == 4) {
            result = derivative4(t);
        } else if (order == 5) {
            result = derivative5(t);
        } else {
            std::string error_msg = std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has 1st, 2nd, 3rd, 4th "
                "order derivatives: order = {0:d}.",
                order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        return result;
    }

    [[using gnu: pure, flatten, leaf]]
    T integral(U lower_bound, U upper_bound) const noexcept {
        constexpr std::array<U, 2> kFactors{-(1.0 / 24.0), 1.0 / 240.0};
        int sign{1};
        if (lower_bound > upper_bound) {
            spdlog::warn(
                "Invalid argument issue detected! the lower_bound of integral should always be "
                "less than upper_bound: lower_bound = {0:.6f} while upper_bound = {1:.6f}.",
                lower_bound, upper_bound
            );
            std::swap(lower_bound, upper_bound);
            sign = -1;
        }
        const std::size_t size{ts_.size()};
        const std::size_t istart =
            std::lower_bound(ts_.cbegin(), ts_.cend(), lower_bound) - ts_.cbegin();
        const std::size_t iend =
            std::lower_bound(ts_.cbegin(), ts_.cend(), upper_bound) - ts_.cbegin();
        U h, h3;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (operator()(lower_bound) + operator()(upper_bound)) * h * 0.5 * sign;
        }
        T result{0.0};
        h = ts_[istart] - lower_bound;
        h3 = h * h * h;
        result += (operator()(lower_bound) + ys_[istart]) * h * 0.5 +
                  (derivative2(lower_bound) + ddys_[istart]) * h3 * kFactors[0] +
                  (derivative4(lower_bound) + d4ys_[istart]) * h3 * h * h * kFactors[1];
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            h = ts_[i + 1] - ts_[i];
            h3 = h * h * h;
            result += (ys_[i] + ys_[i + 1]) * h * 0.5 +
                      (ddys_[i] + ddys_[i + 1]) * h3 * kFactors[0] +
                      (d4ys_[i] + d4ys_[i + 1]) * h3 * h * h * kFactors[1];
        }
        h = upper_bound - ts_[iend - 1];
        h3 = h * h * h;
        result += (ys_[iend - 1] + operator()(upper_bound)) * h * 0.5 +
                  (ddys_[iend - 1] + derivative2(upper_bound)) * h3 * kFactors[0] +
                  (d4ys_[iend - 1] + derivative4(upper_bound)) * h3 * h * h * kFactors[1];
        result *= sign;
        return result;
    }

    [[using gnu: pure, always_inline]]
    U minT() const noexcept {
        return ts_.front();
    }

    [[using gnu: pure, always_inline]]
    U maxT() const noexcept {
        return ts_.back();
    }

    template <typename _T = T, std::enable_if_t<std::is_floating_point_v<_T>, bool> = true>
    [[using gnu: pure, flatten, leaf]]
    T minY() const noexcept {
        std::size_t pos = std::min_element(ys_.cbegin(), ys_.cend()) - ys_.cbegin();
        U h, ratio;
        T derivative, derivative2;
        if (pos == 0) {
            h = ts_[1] - ts_[0];
            derivative = quinerpd(ys_[0], ys_[1], ddys_[0], ddys_[1], ys_[0], ys_[1], 0.0, h);
            if (derivative > 0.0) {
                return ys_[0];
            } else {
                pos += 1;
            }
        } else if (pos == ts_.size() - 1) {
            h = ts_[pos] - ts_[pos - 1];
            derivative = quinerpd(
                ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], 1.0,
                h
            );
            if (derivative < 0.0) {
                return ys_[pos];
            }
        } else {
            h = ts_[pos + 1] - ts_[pos];
            derivative = quinerpd(
                ys_[pos], ys_[pos + 1], ddys_[pos], ddys_[pos + 1], d4ys_[pos], d4ys_[pos + 1], 0.0,
                h
            );
            if (derivative < 0.0) {
                pos += 1;
            }
        }
        h = ts_[pos] - ts_[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter; --num_iter) {
            derivative = quinerpd(
                ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos],
                ratio, h
            );
            derivative2 = cuberp(ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h);
            ratio = derivative / (derivative2 * h);
        }
        return quinerp(
            ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h
        );
    }

    template <typename _T = T, std::enable_if_t<std::is_floating_point_v<_T>, bool> = true>
    [[using gnu: pure, flatten, leaf]]
    T maxY() const noexcept {
        std::size_t pos = std::max_element(ys_.cbegin(), ys_.cend()) - ys_.cbegin();
        U h, ratio;
        T derivative, derivative2;
        if (pos == 0) {
            h = ts_[1] - ts_[0];
            derivative = quinerpd(ys_[0], ys_[1], ddys_[0], ddys_[1], ys_[0], ys_[1], 0.0, h);
            if (derivative < 0.0) {
                return ys_[0];
            } else {
                pos += 1;
            }
        } else if (pos == ts_.size() - 1) {
            h = ts_[pos] - ts_[pos - 1];
            derivative = quinerpd(
                ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], 1.0,
                h
            );
            if (derivative > 0.0) {
                return ys_[pos];
            }
        } else {
            h = ts_[pos + 1] - ts_[pos];
            derivative = quinerpd(
                ys_[pos], ys_[pos + 1], ddys_[pos], ddys_[pos + 1], d4ys_[pos], d4ys_[pos + 1], 0.0,
                h
            );
            if (derivative > 0.0) {
                pos += 1;
            }
        }
        h = ts_[pos] - ts_[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter; --num_iter) {
            derivative = quinerpd(
                ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos],
                ratio, h
            );
            derivative2 = cuberp(ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h);
            ratio = derivative / (derivative2 * h);
        }
        return quinerp(
            ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    const std::vector<U>& ts() const noexcept {
        return ts_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& ys() const noexcept {
        return ys_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& ddys() const noexcept {
        return ddys_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& d4ys() const noexcept {
        return d4ys_;
    }

  private:
    struct [[nodiscard]] OutriggerMatrix final {
        [[using gnu: pure, always_inline]] [[nodiscard]]
        std::vector<T> gaussSeidel(const std::vector<T>& b) const noexcept {
            const std::size_t mat_size{a_diag.size()};
            const std::size_t half_mat_size{mat_size / 2};
            std::vector<T> x(mat_size, 0.0);
            for (std::size_t num_iter{40}; num_iter; --num_iter) {
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

    [[using gnu: flatten, leaf]] explicit PiecewiseQuinticFunction1(
        std::vector<U> ts, std::vector<T> ys, std::vector<T> ddys, std::vector<T> d4ys
    ) {
        if (ts.size() < 2 || ys.size() < 2) {
            std::string error_msg = std::format(
                "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() = "
                "{0:d} while ys.size() = {1:d}",
                ts.size(), ys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        } else if (ts.size() != ys.size()) {
            std::string error_msg = std::format(
                "Invalid arguments detected! ts, ys must share the same size: ts.size() = {0:d} "
                "while ys.size() = {1:d}",
                ts.size(), ys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        } else if (!std::is_sorted(ts.cbegin(), ts.cend())) {
            std::string error_msg =
                std::format("Invalid arguments detected! ts has to be a sorted array!");
            throw std::invalid_argument(std::move(error_msg));
        } else if (hasDuplicates(ts.cbegin(), ts.cend(), kDuplicateCriterion)) {
            std::string error_msg =
                std::format("Invalid arguments detected! ts can not have duplicated elements!");
            throw std::invalid_argument(std::move(error_msg));
        }
        if (ddys.size() != ts.size()) {
            std::string error_msg = std::format(
                "Invalid arguments detected! ts, ys, ddys must share the same size: ts.size() = "
                "{0:d} while ddys.size() = {1:d}",
                ts.size(), ddys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (d4ys.size() != ts.size()) {
            std::string error_msg = std::format(
                "Invalid arguments detected! ts, ys, ddys must share the same size: ts.size() = "
                "{0:d} while d4ys.size() = {1:d}",
                ts.size(), d4ys.size()
            );
            throw std::invalid_argument(std::move(error_msg));
        }

        ts_ = std::move(ts);
        ys_ = std::move(ys);
        ddys_ = std::move(ddys);
        d4ys_ = std::move(d4ys);
    }

    [[using gnu: pure, always_inline]]
    T derivative2(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return T{0.0};
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return T{0.0};
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio = (t - ts_[pos - 1]) / h;
        return cuberp(ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    T derivative3(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return T{0.0};
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return T{0.0};
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio{(t - ts_[pos - 1]) / h};
        return cuberpd(ddys_[pos - 1], ddys_[pos], d4ys_[pos - 1], d4ys_[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    T derivative4(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return T{0.0};
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return T{0.0};
        }
        const U ratio{(t - ts_[pos - 1]) / (ts_[pos] - ts_[pos - 1])};
        return lerp(d4ys_[pos - 1], d4ys_[pos], ratio);
    }

    [[using gnu: pure, always_inline]]
    T derivative5(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return T{0.0};
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return T{0.0};
        }
        return (d4ys_[pos] - d4ys_[pos - 1]) / (ts_[pos] - ts_[pos - 1]);
    }

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& ts_;
        ar& ys_;
        ar& ddys_;
        ar& d4ys_;
        return;
    }

    std::vector<U> ts_;
    std::vector<T> ys_;
    std::vector<T> ddys_;
    std::vector<T> d4ys_;
};

template <typename T, typename U>
struct isFunction<PiecewiseQuinticFunction1<T, U>> final {
    static constexpr bool value = true;
};

using PiecewiseQuinticFunction1f = PiecewiseQuinticFunction1<float>;
using PiecewiseQuinticFunction1d = PiecewiseQuinticFunction1<double>;

} // namespace math
} // namespace tiny_pnc
