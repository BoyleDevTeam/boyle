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

PRO_DEF_MEM_DISPATCH(m_eval, eval);
PRO_DEF_MEM_DISPATCH(m_tangent, tangent);
PRO_DEF_MEM_DISPATCH(m_orthonormal, orthonormal);
PRO_DEF_MEM_DISPATCH(m_inverse, inverse);
PRO_DEF_MEM_DISPATCH(m_minS, minS);
PRO_DEF_MEM_DISPATCH(m_maxS, maxS);
PRO_DEF_MEM_DISPATCH(m_front, front);
PRO_DEF_MEM_DISPATCH(m_back, back);
PRO_DEF_MEM_DISPATCH(m_arcLengths, arcLengths);
PRO_DEF_MEM_DISPATCH(m_anchorPoints, anchorPoints);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <std::floating_point T = double>
class [[nodiscard]] Curve2Facade final : public pro::facade_builder
    ::template add_convention<m_eval, auto(T) noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<m_tangent, auto(T) noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<m_orthonormal, auto(T) noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<m_inverse, auto(::boyle::math::Vec2<T>) noexcept -> ::boyle::math::SlDuplet<T>>
    ::template add_convention<m_minS, auto() noexcept -> T>
    ::template add_convention<m_maxS, auto() noexcept -> T>
    ::template add_convention<m_front, auto() noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<m_back, auto() noexcept -> ::boyle::math::Vec2<T>>
    ::template add_convention<m_arcLengths, auto() noexcept -> const std::vector<T>&>
    ::template add_convention<m_anchorPoints, auto() noexcept -> const std::vector<::boyle::math::Vec2<T>>&>
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
