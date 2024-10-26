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
#include <type_traits>
#include <utility>
#include <vector>

#include "proxy.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/duplet.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_tangent, tangent);
PRO_DEF_MEM_DISPATCH(c_normal, normal);
PRO_DEF_MEM_DISPATCH(c_curvature, curvature);
PRO_DEF_MEM_DISPATCH(c_inverse, inverse);
PRO_DEF_MEM_DISPATCH(c_minS, minS);
PRO_DEF_MEM_DISPATCH(c_maxS, maxS);
PRO_DEF_MEM_DISPATCH(c_front, front);
PRO_DEF_MEM_DISPATCH(c_back, back);
PRO_DEF_MEM_DISPATCH(c_arcLengths, arcLengths);
PRO_DEF_MEM_DISPATCH(c_anchorPoints, anchorPoints);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <VecArithmetic T, std::floating_point U = typename T::value_type>
class [[nodiscard]] Curve2Facade final : public pro::facade_builder
    ::template add_convention<pro::operator_dispatch<"()">, auto(U) const noexcept -> T>
    ::template add_convention<c_eval, auto(U) const noexcept -> T>
    ::template add_convention<c_tangent, auto(U) const noexcept -> T>
    ::template add_convention<c_normal, auto(U) const noexcept -> T>
    ::template add_convention<c_curvature, auto(U) const noexcept -> U>
    ::template add_convention<c_inverse, auto(T) const noexcept -> ::boyle::math::SlDuplet<U>>
    ::template add_convention<c_minS, auto() const noexcept -> U>
    ::template add_convention<c_maxS, auto() const noexcept -> U>
    ::template add_convention<c_front, auto() const noexcept -> T>
    ::template add_convention<c_back, auto() const noexcept -> T>
    ::template add_convention<c_arcLengths, auto() const noexcept -> const std::vector<U>&>
    ::template add_convention<c_anchorPoints, auto() const noexcept -> const std::vector<T>&>
    ::build {};

// clang-format on

} // namespace detail

template <VecArithmetic T, std::floating_point U = typename T::value_type>
using Curve2Proxy = pro::proxy<detail::Curve2Facade<T, U>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeCurve2Proxy(T&& curve2) -> Curve2Proxy<
                                            typename std::remove_reference_t<T>::value_type,
                                            typename std::remove_reference_t<T>::param_type> {
    return pro::make_proxy<detail::Curve2Facade<
        typename std::remove_reference_t<T>::value_type,
        typename std::remove_reference_t<T>::param_type>>(std::forward<T>(curve2));
}

} // namespace boyle::math
