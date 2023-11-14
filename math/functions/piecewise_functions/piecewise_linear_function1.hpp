/**
 * @file piecewise_linear_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-10
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <format>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "spdlog/spdlog.h"

#include "common/utils/macros.hpp"
#include "math/functions/function1.hpp"
#include "math/type_traits.hpp"
#include "math/utils.hpp"

namespace tiny_pnc {
namespace math {

template <typename T, typename U = T>
class [[nodiscard]] PiecewiseLinearFunction1 final
    : public Function1<PiecewiseLinearFunction1<T, U>, T, U> {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    static constexpr U kDuplicateCriterion{kEpsilon};

    ENABLE_IMPLICIT_CONSTRUCTORS(PiecewiseLinearFunction1);

    [[using gnu: flatten,
      leaf]] explicit PiecewiseLinearFunction1(std::vector<U> ts, std::vector<T> ys) {
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
    }

    ~PiecewiseLinearFunction1() noexcept override = default;

    [[using gnu: pure, always_inline]] T operator()(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return lerp(ys_[0], ys_[1], (t - ts_[0]) / (ts_[1] - ts_[0]));
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return lerp(
                ys_[pos - 1], ys_[pos - 2], (t - ts_[pos - 1]) / (ts_[pos - 2] - ts_[pos - 1])
            );
        }
        return lerp(ys_[pos - 1], ys_[pos], (t - ts_[pos - 1]) / (ts_[pos] - ts_[pos - 1]));
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t) const noexcept {
        const std::size_t pos = std::upper_bound(ts_.cbegin(), ts_.cend(), t) - ts_.cbegin();
        if (pos == 0) {
            spdlog::warn(
                "Out of range issue detected! t should be greater than minT(): t = {0:.6f} while "
                "minT() = {1:.6f}. Use extra interpolating!",
                t, ts_.front()
            );
            return (ys_[1] - ys_[0]) / (ts_[1] - ts_[0]);
        } else if (pos == ts_.size()) {
            spdlog::warn(
                "Out of range issue detected! t should be less than maxT(): t = {0:.6f} while "
                "maxT() = {1:.6f}. Use extra interpolating",
                t, ts_.back()
            );
            return (ys_[pos - 1] - ys_[pos - 2]) / (ts_[pos - 1] - ts_[pos - 2]);
        }
        return (ys_[pos] - ys_[pos - 1]) / (ts_[pos] - ts_[pos - 1]);
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t, unsigned int order) const {
        if (order != 1) {
            std::string error_msg = std::format(
                "Invalid argument error! The PiecewiseLinearFunction only has first order "
                "derivative: order = {0:d}.",
                order
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        return derivative(t);
    }

    [[using gnu: pure]]
    T integral(U lower_bound, U upper_bound) const noexcept {
        int sign{1};
        if (lower_bound > upper_bound) {
            spdlog::warn(
                "Invalid argument issue detected! The lower_bound of integral should always be "
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
        if (istart == size || iend == 0 || istart == iend) {
            return (operator()(lower_bound) + operator()(upper_bound)) *
                   (upper_bound - lower_bound) * 0.5 * sign;
        }
        T result{0.0};
        result += (operator()(lower_bound) + ys_[istart]) * (ts_[istart] - lower_bound) * 0.5;
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            result += (ys_[i] + ys_[i + 1]) * (ts_[i + 1] - ts_[i]) * 0.5;
        }
        result += (ys_[iend - 1] + operator()(upper_bound)) * (upper_bound - ts_[iend - 1]) * 0.5;
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
    [[using gnu: pure, always_inline]]
    T minY() const noexcept {
        return *std::min_element(ys_.cbegin(), ys_.cend());
    }

    template <typename _T = T, std::enable_if_t<std::is_floating_point_v<_T>, bool> = true>
    [[using gnu: pure, always_inline]]
    T maxY() const noexcept {
        return *std::max_element(ys_.cbegin(), ys_.cend());
    }

    [[using gnu: pure, always_inline]]
    const std::vector<U>& ts() const noexcept {
        return ts_;
    }

    [[using gnu: pure, always_inline]]
    const std::vector<T>& ys() const noexcept {
        return ys_;
    }

  private:
    std::vector<U> ts_;
    std::vector<T> ys_;

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& ts_;
        ar& ys_;
        return;
    }
};

template <typename T, typename U>
struct isFunction<PiecewiseLinearFunction1<T, U>> final {
    static constexpr bool value = true;
};

using PiecewiseLinearFunction1f = PiecewiseLinearFunction1<float>;
using PiecewiseLinearFunction1d = PiecewiseLinearFunction1<double>;

} // namespace math
} // namespace tiny_pnc
