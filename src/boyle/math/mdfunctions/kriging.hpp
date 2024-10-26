/**
 * @file kriging.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-16
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <functional>
#include <span>

namespace boyle::math {

template <std::floating_point T, typename U = std::span<const T>>
class [[nodiscard]] Kriging final {
  public:
    using value_type = T;
    using param_type = U;

    Kriging() noexcept = default;
    Kriging(const Kriging& other) noexcept = default;
    auto operator=(const Kriging& other) noexcept -> Kriging& = default;
    Kriging(Kriging&& other) noexcept = default;
    auto operator=(Kriging&& other) noexcept -> Kriging& = default;
    ~Kriging() noexcept = default;

    [[nodiscard]]
    auto num_dimensions() const noexcept -> std::size_t;
    [[nodiscard]]
    auto num_points() const noexcept -> std::size_t;

    auto eval(param_type x) const noexcept -> value_type;

    auto gradient(param_type x) const noexcept -> std::vector<value_type>;

    auto gradient(param_type x, std::size_t idx) const noexcept -> value_type;

    [[nodiscard]]
    auto has_extrema(param_type x) const noexcept -> bool;

  private:
    std::function<auto(param_type) noexcept -> value_type> m_variogram;
};
} // namespace boyle::math
