/**
 * @file piecewise_cubic_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-15
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
#include "math/type_traits.hpp"
#include "math/utils.hpp"

namespace tiny_pnc {
namespace math {

template <typename T, typename U = T>
class [[nodiscard]] PiecewiseCubicFunction1 final
    : public Function1<PiecewiseCubicFunction1<T, U>, T, U> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    static constexpr U kDuplicateCriterion{kEpsilon};

    struct BoundaryMode final {
        unsigned int order;
        T derivative;
    };

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseCubicFunction1);

    [[using gnu: always_inline]] explicit PiecewiseCubicFunction1(
        std::vector<U> ts, std::vector<T> ys
    )
        : PiecewiseCubicFunction1<T, U>(
              std::move(ts), std::move(ys), BoundaryMode{2, 0.0}, BoundaryMode{2, 0.0}
          ) {}

    [[using gnu: flatten, leaf]] explicit PiecewiseCubicFunction1(
        std::vector<U> ts, std::vector<T> ys, BoundaryMode b0, BoundaryMode bf
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

        ts_ = std::move(ts);
        ys_ = std::move(ys);

        if (b0.order == 0 || b0.order > 2) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "b0.order = {0:d}.",
                b0.order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (bf.order == 0 || bf.order > 2) {
            std::string error_msg = std::format(
                "Invalid argument detected! The derivative order of b0 can only be 1, 2, 3: "
                "bf.order = {0:d}.",
                bf.order
            );
            throw std::invalid_argument(std::move(error_msg));
        }

        const std::size_t size{ts_.size()};
        std::vector<U> hs(size - 1);
        std::vector<T> ds(size - 1);
        std::vector<U> a_diag(size);
        std::vector<T> b(size);

        hs[0] = ts_[1] - ts_[0];
        ds[0] = (ys_[1] - ys_[0]) / hs[0];
        for (std::size_t i{1}; i < size - 1; i++) {
            hs[i] = ts_[i + 1] - ts_[i];
            ds[i] = (ys_[i + 1] - ys_[i]) / hs[i];
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

        ddys_ = A.luDcmp(b);
    }

    ~PiecewiseCubicFunction1() noexcept override = default;

    [[using gnu: pure, always_inline]] T operator()(U t) const noexcept {
        constexpr std::array<U, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
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
                   (ddys_[0] * kFactors[0] + ddys_[1] * kFactors[1]) * (t - ts_[0]) * h;
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
                       (t - ts_[pos - 1]) * h;
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio{(t - ts_[pos - 1]) / (ts_[pos] - ts_[pos - 1])};
        return cuberp(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t) const noexcept {
        constexpr std::array<U, 2> kFactors{-(1.0 / 3.0), -(1.0 / 6.0)};
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            const U h{ts_[1] - ts_[0]};
            return (ys_[1] - ys_[0]) / h + (ddys_[0] * kFactors[0] + ddys_[1] * kFactors[1]) * h;
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            const U h{ts_[pos - 2] - ts_[pos - 1]};
            return (ys_[pos - 2] - ys_[pos - 1]) / h +
                   (ddys_[pos - 1] * kFactors[0] + ddys_[pos - 2] * kFactors[1]) * h;
        }
        const U h{ts_[pos] - ts_[pos - 1]};
        const U ratio{(t - ts_[pos - 1]) / h};
        return cuberpd(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
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
        } else {
            std::string error_msg = std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has 1, 2, 3 order "
                "derivatives: order = {0:d}.",
                order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        return result;
    }

    [[using gnu: pure]]
    T integral(U lower_bound, U upper_bound) const noexcept {
        constexpr U kFactor{-(1.0 / 24.0)};
        int sign = 1;
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
        U h;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (operator()(lower_bound) + operator()(upper_bound)) * h * 0.5 * sign;
        }
        T result{0.0};
        h = ts_[istart] - lower_bound;
        result += (operator()(lower_bound) + ys_[istart]) * h * 0.5 +
                  (derivative2(lower_bound) + ddys_[istart]) * h * h * h * kFactor;
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            h = ts_[i + 1] - ts_[i];
            result +=
                (ys_[i] + ys_[i + 1]) * h * 0.5 + (ddys_[i] + ddys_[i + 1]) * h * h * h * kFactor;
        }
        h = upper_bound - ts_[iend - 1];
        result += (ys_[iend - 1] + operator()(upper_bound)) * h * 0.5 +
                  (ddys_[iend - 1] + derivative2(upper_bound)) * h * h * h * kFactor;
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
            derivative = cuberpd(ys_[0], ys_[1], ddys_[0], ddys_[1], 0.0, h);
            if (derivative > 0.0) {
                return ys_[0];
            } else {
                pos += 1;
            }
        } else if (pos == ts_.size() - 1) {
            h = ts_[pos] - ts_[pos - 1];
            derivative = cuberpd(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], 1.0, h);
            if (derivative < 0.0) {
                return ys_[pos];
            }
        } else {
            h = ts_[pos + 1] - ts_[pos];
            derivative = cuberpd(ys_[pos], ys_[pos + 1], ddys_[pos], ddys_[pos + 1], 0.0, h);
            if (derivative < 0.0) {
                pos += 1;
            }
        }
        h = ts_[pos] - ts_[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter; --num_iter) {
            derivative = cuberpd(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
            derivative2 = lerp(ddys_[pos - 1], ddys_[pos], ratio);
            ratio -= derivative / (derivative2 * h);
        }
        return cuberp(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
    }

    template <typename _T = T, std::enable_if_t<std::is_floating_point_v<_T>, bool> = true>
    [[using gnu: pure, flatten, leaf]]
    T maxY() const noexcept {
        std::size_t pos = std::max_element(ys_.cbegin(), ys_.cend()) - ys_.cbegin();
        U h, ratio;
        T derivative, derivative2;
        if (pos == 0) {
            h = ts_[1] - ts_[0];
            derivative = cuberpd(ys_[0], ys_[1], ddys_[0], ddys_[1], 0.0, h);
            if (derivative < 0.0) {
                return ys_[0];
            } else {
                pos += 1;
            }
        } else if (pos == ts_.size() - 1) {
            h = ts_[pos] - ts_[pos - 1];
            derivative = cuberpd(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], 1.0, h);
            if (derivative > 0.0) {
                return ys_[pos];
            }
        } else {
            h = ts_[pos + 1] - ts_[pos];
            derivative = cuberpd(ys_[pos], ys_[pos + 1], ddys_[pos], ddys_[pos + 1], 0.0, h);
            if (derivative > 0.0) {
                pos += 1;
            }
        }
        h = ts_[pos] - ts_[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter; --num_iter) {
            derivative = cuberpd(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
            derivative2 = lerp(ddys_[pos - 1], ddys_[pos], ratio);
            ratio -= derivative / (derivative2 * h);
        }
        return cuberp(ys_[pos - 1], ys_[pos], ddys_[pos - 1], ddys_[pos], ratio, h);
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

  private:
    struct [[nodiscard]] TridiagonalMatrix final {
        [[using gnu: pure, always_inline]] [[nodiscard]]
        std::vector<T> luDcmp(const std::vector<T>& b) const noexcept {
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

    [[using gnu: flatten, leaf]] explicit PiecewiseCubicFunction1(
        std::vector<U> ts, std::vector<T> ys, std::vector<T> ddys
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

        ts_ = std::move(ts);
        ys_ = std::move(ys);
        ddys_ = std::move(ddys);
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
        const U ratio = (t - ts_[pos - 1]) / (ts_[pos] - ts_[pos - 1]);
        return lerp(ddys_[pos - 1], ddys_[pos], ratio);
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
        return (ddys_[pos] - ddys_[pos - 1]) / (ts_[pos] - ts_[pos - 1]);
    }

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& ts_;
        ar& ys_;
        ar& ddys_;
        return;
    }

    std::vector<U> ts_;
    std::vector<T> ys_;
    std::vector<T> ddys_;
};

template <typename T, typename U>
struct isFunction<PiecewiseCubicFunction1<T, U>> final {
    static constexpr bool value = true;
};

using PiecewiseCubicFunction1f = PiecewiseCubicFunction1<float>;
using PiecewiseCubicFunction1d = PiecewiseCubicFunction1<double>;

} // namespace math
} // namespace tiny_pnc
