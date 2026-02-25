/**
 * @file curve_proxy.hpp
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
#include <span>
#include <utility>

#include "proxy/proxy.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/curves/sl.hpp"
#include "boyle/math/curves/slv.hpp"
#include "boyle/math/dense/vec2.hpp"
#include "boyle/math/dense/vec3.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_tangent, tangent);
PRO_DEF_MEM_DISPATCH(c_normal, normal);
PRO_DEF_MEM_DISPATCH(c_binormal, binormal);
PRO_DEF_MEM_DISPATCH(c_curvature, curvature);
PRO_DEF_MEM_DISPATCH(c_torsion, torsion);
PRO_DEF_MEM_DISPATCH(c_inverse, inverse);
PRO_DEF_MEM_DISPATCH(c_minS, minS);
PRO_DEF_MEM_DISPATCH(c_maxS, maxS);
PRO_DEF_MEM_DISPATCH(c_front, front);
PRO_DEF_MEM_DISPATCH(c_back, back);
PRO_DEF_MEM_DISPATCH(c_arcLengths, arcLengths);
PRO_DEF_MEM_DISPATCH(c_anchorPoints, anchorPoints);

// NOLINTEND(modernize-use-trailing-return-type)

template <typename T>
class CurveFacade final {
    static_assert(false, "::boyle::math::CurveFacade is not implemented for this type argument!");
};

// clang-format off

template <Vec2Arithmetic T>
    requires std::floating_point<typename T::value_type>
class CurveFacade<T> final : public pro::facade_builder
    ::support_copy<pro::constraint_level::nontrivial>
    ::template add_convention<
        pro::operator_dispatch<"()">,
        auto(typename T::value_type) const noexcept -> T,
        auto(typename T::value_type, typename T::value_type) const noexcept -> T,
        auto(::boyle::math::SlDuplet<typename T::value_type>) const noexcept -> T
    >
    ::template add_convention<
        c_eval,
        auto(typename T::value_type) const noexcept -> T,
        auto(typename T::value_type, typename T::value_type) const noexcept -> T,
        auto(::boyle::math::SlDuplet<typename T::value_type>) const noexcept -> T
    >
    ::template add_convention<c_tangent, auto(typename T::value_type) const noexcept -> T>
    ::template add_convention<c_normal, auto(typename T::value_type) const noexcept -> T>
    ::template add_convention<c_curvature, auto(typename T::value_type) const noexcept -> typename T::value_type>
    ::template add_convention<
        c_inverse,
        auto(T) const noexcept -> ::boyle::math::SlDuplet<typename T::value_type>, 
        auto(T, typename T::value_type, typename T::value_type) const noexcept -> ::boyle::math::SlDuplet<typename T::value_type>
    >
    ::template add_convention<c_minS, auto() const noexcept -> typename T::value_type>
    ::template add_convention<c_maxS, auto() const noexcept -> typename T::value_type>
    ::template add_convention<c_front, auto() const noexcept -> T>
    ::template add_convention<c_back, auto() const noexcept -> T>
    ::template add_convention<c_arcLengths, auto() const noexcept -> std::span<const typename T::value_type>>
    ::template add_convention<c_anchorPoints, auto() const noexcept -> std::span<const T>>
    ::build {};

template <Vec3Arithmetic T>
    requires std::floating_point<typename T::value_type>
class CurveFacade<T> final : public pro::facade_builder
    ::support_copy<pro::constraint_level::nontrivial>
    ::template add_convention<
        pro::operator_dispatch<"()">,
        auto(typename T::value_type) const noexcept -> T,
        auto(typename T::value_type, typename T::value_type, typename T::value_type) const noexcept -> T,
        auto(::boyle::math::SlvTriplet<typename T::value_type>) const noexcept -> T
    >
    ::template add_convention<
        c_eval,
        auto(typename T::value_type) const noexcept -> T,
        auto(typename T::value_type, typename T::value_type, typename T::value_type) const noexcept -> T,
        auto(::boyle::math::SlvTriplet<typename T::value_type>) const noexcept -> T
    >
    ::template add_convention<c_tangent, auto(typename T::value_type) const noexcept -> T>
    ::template add_convention<c_normal, auto(typename T::value_type) const noexcept -> T>
    ::template add_convention<c_binormal, auto(typename T::value_type) const noexcept -> T>
    ::template add_convention<c_curvature, auto(typename T::value_type) const noexcept -> typename T::value_type>
    ::template add_convention<c_torsion, auto(typename T::value_type) const noexcept -> typename T::value_type>
    ::template add_convention<
        c_inverse,
        auto(T) const noexcept -> ::boyle::math::SlvTriplet<typename T::value_type>,
        auto(T, typename T::value_type, typename T::value_type) const noexcept -> ::boyle::math::SlvTriplet<typename T::value_type>
    >
    ::template add_convention<c_minS, auto() const noexcept -> typename T::value_type>
    ::template add_convention<c_maxS, auto() const noexcept -> typename T::value_type>
    ::template add_convention<c_front, auto() const noexcept -> T>
    ::template add_convention<c_back, auto() const noexcept -> T>
    ::template add_convention<c_arcLengths, auto() const noexcept -> std::span<const typename T::value_type>>
    ::template add_convention<c_anchorPoints, auto() const noexcept -> std::span<const T>>
    ::build {};

// clang-format on

} // namespace detail

template <VecArithmetic T>
using CurveProxy = pro::proxy<detail::CurveFacade<T>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeCurveProxy(T&& curve) -> CurveProxy<typename std::remove_cvref_t<T>::value_type> {
    return pro::make_proxy<detail::CurveFacade<typename std::remove_cvref_t<T>::value_type>>(
        std::forward<T>(curve)
    );
}

} // namespace boyle::math
