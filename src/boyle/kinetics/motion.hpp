/**
 * @file motion.hpp
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

#include <array>
#include <memory_resource>
#include <span>

#include "boost/serialization/access.hpp"

#include "boyle/math/dense/detail/dense_degenerate_trait.hpp"
#include "boyle/math/functions/boundary_mode.hpp"
#include "boyle/math/functions/piecewise_quintic_function.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class Motion final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = ::boyle::math::detail::DenseDegenerateTraitT<value_type>;
    using size_type = std::size_t;
    using allocator_type = std::pmr::polymorphic_allocator<value_type>;

    Motion() noexcept = default;
    Motion(const Motion& other) noexcept = default;
    auto operator=(const Motion& other) noexcept -> Motion& = default;
    Motion(Motion&& other) noexcept = default;
    auto operator=(Motion&& other) noexcept -> Motion& = default;
    ~Motion() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_s_of_t.get_allocator();
    }

    template <
        std::ranges::input_range R0 = std::initializer_list<param_type>,
        std::ranges::input_range R1 = std::initializer_list<value_type>>
    [[using gnu: always_inline]]
    explicit Motion(R0&& ts, R1&& ss, const allocator_type& alloc = {})
        requires std::same_as<std::ranges::range_value_t<R0>, param_type> &&
                 std::same_as<std::ranges::range_value_t<R1>, value_type>
        : Motion(
              ts, ss,
              std::array<::boyle::math::BoundaryMode<value_type>, 2>{
                  ::boyle::math::BoundaryMode<value_type>{2, value_type{0.0}},
                  ::boyle::math::BoundaryMode<value_type>{4, value_type{0.0}}
              },
              std::array<::boyle::math::BoundaryMode<value_type>, 2>{
                  ::boyle::math::BoundaryMode<value_type>{2, value_type{0.0}},
                  ::boyle::math::BoundaryMode<value_type>{4, value_type{0.0}}
              },
              alloc
          ) {}

    template <
        std::ranges::input_range R0 = std::initializer_list<param_type>,
        std::ranges::input_range R1 = std::initializer_list<value_type>>
    [[using gnu: always_inline]]
    explicit Motion(
        R0&& ts, R1&& ss, std::array<::boyle::math::BoundaryMode<value_type>, 2> b0,
        std::array<::boyle::math::BoundaryMode<value_type>, 2> bf, const allocator_type& alloc = {}
    )
        requires std::same_as<std::ranges::range_value_t<R0>, param_type> &&
                 std::same_as<std::ranges::range_value_t<R1>, value_type>
        : m_s_of_t(ts, ss, b0, bf, alloc) {}

    [[using gnu: pure, always_inline]]
    auto s(param_type t) const noexcept -> value_type {
        return m_s_of_t.eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto velocity(param_type t) const noexcept -> value_type {
        return m_s_of_t.derivative(t);
    }

    [[using gnu: pure, always_inline]]
    auto accel(param_type t) const noexcept -> value_type {
        return m_s_of_t.derivative(t, 2);
    }

    [[using gnu: pure, always_inline]]
    auto jerk(param_type t) const noexcept -> value_type {
        return m_s_of_t.derivative(t, 3);
    }

    [[using gnu: pure, always_inline]]
    auto snap(param_type t) const noexcept -> value_type {
        return m_s_of_t.derivative(t, 4);
    }

    [[using gnu: pure, always_inline]]
    auto minT() const noexcept -> param_type {
        return m_s_of_t.minT();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> param_type {
        return m_s_of_t.maxT();
    }

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> std::span<const param_type> {
        return m_s_of_t.ts();
    }

    [[using gnu: pure, always_inline]]
    auto ss() const noexcept -> std::span<const value_type> {
        return m_s_of_t.ys();
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_s_of_t;
        return;
    }

    ::boyle::math::PiecewiseQuinticFunction<value_type, allocator_type> m_s_of_t;
};

using Motion1s = Motion<float>;

using Motion1d = Motion<double>;

} // namespace boyle::kinetics
