/**
 * @file motion1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <type_traits>
#include <vector>

#include "boost/serialization/access.hpp"

#include "common/utils/macros.hpp"
#include "math/functions/piecewise_functions/piecewise_quintic_function1.hpp"

namespace tiny_pnc {
namespace kinetics {

template <typename T>
class [[nodiscard]] Motion1 final {
    static_assert(std::is_floating_point_v<T>, "The loaded type must be a floating-point type.");
    friend class boost::serialization::access;

  public:
    using BoundaryMode = typename tiny_pnc::math::PiecewiseQuinticFunction1<T>::BoundaryMode;

    [[using gnu: always_inline]] explicit Motion1(std::vector<T> ts, std::vector<T> ss)
        : Motion1(
              std::move(ts), std::move(ss),
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}},
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}}
          ) {}

    [[using gnu: always_inline]] explicit Motion1(
        std::vector<T> ts, std::vector<T> ss, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf
    )
        : s_of_t_(std::move(ts), std::move(ss), b0, bf) {}

    ENABLE_IMPLICIT_CONSTRUCTORS(Motion1);

    ~Motion1() noexcept = default;

    [[using gnu: pure, always_inline]]
    T s(T t) const noexcept {
        return s_of_t_(t);
    }

    [[using gnu: pure, always_inline]]
    T velocity(T t) const noexcept {
        return s_of_t_.derivative(t);
    }

    [[using gnu: pure, always_inline]]
    T accel(T t) const noexcept {
        return s_of_t_.derivative(t, 2);
    }

    [[using gnu: pure, always_inline]]
    T jerk(T t) const noexcept {
        return s_of_t_.derivative(t, 3);
    }

    [[using gnu: pure, always_inline]]
    T snap(T t) const noexcept {
        return s_of_t_.derivative(t, 4);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> s(const std::vector<T>& ts) const noexcept {
        return s_of_t_(ts);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> velocity(const std::vector<T>& ts) const noexcept {
        return s_of_t_.derivative(ts);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> accel(const std::vector<T>& ts) const noexcept {
        return s_of_t_.derivative(ts, 2);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> jerk(const std::vector<T>& ts) const noexcept {
        return s_of_t_.derivative(ts, 3);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> snap(const std::vector<T>& ts) const noexcept {
        return s_of_t_.derivative(ts, 4);
    }

  private:
    Motion1(math::PiecewiseQuinticFunction1<T> s_of_t) noexcept : s_of_t_(s_of_t) {}

    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& s_of_t_;
        return;
    }

    tiny_pnc::math::PiecewiseQuinticFunction1<T> s_of_t_;
};

using Motion1f = Motion1<float>;
using Motion1d = Motion1<double>;

} // namespace kinetics
} // namespace tiny_pnc
