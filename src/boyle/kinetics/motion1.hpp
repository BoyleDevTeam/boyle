/**
 * @file motion1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "boyle/math/functions/piecewise_quintic_function1.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class [[nodiscard]] Motion1 final {
    friend class boost::serialization::access;

  public:
    using BoundaryMode = typename ::boyle::math::PiecewiseQuinticFunction1<T>::BoundaryMode;
    ~Motion1() noexcept = default;

    Motion1() noexcept = default;
    Motion1(const Motion1& other) noexcept = default;
    auto operator=(const Motion1& other) noexcept -> Motion1& = default;
    Motion1(Motion1&& other) noexcept = default;
    auto operator=(Motion1&& other) noexcept -> Motion1& = default;

    [[using gnu: always_inline]]
    explicit Motion1(std::vector<T> ts, std::vector<T> ss)
        : Motion1{
              std::move(ts), std::move(ss),
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}},
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}}
          } {}

    [[using gnu: always_inline]]
    explicit Motion1(
        std::vector<T> ts, std::vector<T> ss, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf
    )
        : m_s_of_t{std::move(ts), std::move(ss), b0, bf} {}

    [[using gnu: pure, always_inline]]
    auto s(T t) const noexcept -> T {
        return m_s_of_t.eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto velocity(T t) const noexcept -> T {
        return m_s_of_t.derivative(t);
    }

    [[using gnu: pure, always_inline]]
    auto accel(T t) const noexcept -> T {
        return m_s_of_t.derivative(t, 2);
    }

    [[using gnu: pure, always_inline]]
    auto jerk(T t) const noexcept -> T {
        return m_s_of_t.derivative(t, 3);
    }

    [[using gnu: pure, always_inline]]
    auto snap(T t) const noexcept -> T {
        return m_s_of_t.derivative(t, 4);
    }

    [[using gnu: pure, always_inline]]
    auto minT() const noexcept -> T {
        return m_s_of_t.minT();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> T {
        return m_s_of_t.maxT();
    }

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> std::span<const T> {
        return m_s_of_t.ts();
    }

    [[using gnu: pure, always_inline]]
    auto ss() const noexcept -> std::span<const T> {
        return m_s_of_t.ys();
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_s_of_t;
        return;
    }

    ::boyle::math::PiecewiseQuinticFunction1<T> m_s_of_t{};
};

using Motion1f = Motion1<float>;
using Motion1d = Motion1<double>;

} // namespace boyle::kinetics
