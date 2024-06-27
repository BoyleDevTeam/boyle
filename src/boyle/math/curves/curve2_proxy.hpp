/**
 * @file curve2_proxy.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-24
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <utility>
#include <vector>

#include "proxy.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/utils.hpp"
#include "boyle/math/vec2.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_tangent, tangent);
PRO_DEF_MEM_DISPATCH(c_orthonormal, orthonormal);
PRO_DEF_MEM_DISPATCH(c_inverse, inverse);
PRO_DEF_MEM_DISPATCH(c_minS, minS);
PRO_DEF_MEM_DISPATCH(c_maxS, maxS);
PRO_DEF_MEM_DISPATCH(c_front, front);
PRO_DEF_MEM_DISPATCH(c_back, back);
PRO_DEF_MEM_DISPATCH(c_arcLengths, arcLengths);
PRO_DEF_MEM_DISPATCH(c_anchorPoints, anchorPoints);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <std::floating_point T = double>
class [[nodiscard]] Curve2Facade final : public pro::facade_builder
    ::template add_convention<c_eval, auto(T) const noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<c_tangent, auto(T) const noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<c_orthonormal, auto(T) const noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<c_inverse, auto(::boyle::math::Vec2<T>) const noexcept -> ::boyle::math::SlDuplet<T>>
    ::template add_convention<c_minS, auto() const noexcept -> T>
    ::template add_convention<c_maxS, auto() const noexcept -> T>
    ::template add_convention<c_front, auto() const noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<c_back, auto() const noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<c_arcLengths, auto() const noexcept -> const std::vector<T>&>
    ::template add_convention<c_anchorPoints, auto() const noexcept -> const std::vector<::boyle::math::Vec2<T>>&>
    ::build {};

// clang-format on

} // namespace detail

template <std::floating_point T>
using Curve2Proxy = pro::proxy<detail::Curve2Facade<T>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeCurve2Proxy(T curve2) -> Curve2Proxy<typename T::param_type> {
    return pro::make_proxy<detail::Curve2Facade<typename T::param_type>>(std::move(curve2));
}

} // namespace boyle::math
