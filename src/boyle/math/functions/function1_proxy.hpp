/**
 * @file function1_proxy.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-09
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <utility>

#include "proxy.h"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(m_eval, eval);
PRO_DEF_MEM_DISPATCH(m_derivative, derivative);
PRO_DEF_MEM_DISPATCH(m_integral, integral);
PRO_DEF_MEM_DISPATCH(m_minT, minT);
PRO_DEF_MEM_DISPATCH(m_maxT, maxT);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <::boyle::math::GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] Function1Facader final : public pro::facade_builder
    ::template add_convention<m_eval, auto(U) noexcept -> T>
    ::template add_convention<m_derivative, auto(U) noexcept -> T>
    ::template add_convention<m_integral, auto(U, U) noexcept -> T>
    ::template add_convention<m_minT, auto() noexcept -> U>
    ::template add_convention<m_maxT, auto() noexcept -> U>
    ::build {};

// clang-format on

} // namespace detail

template <GeneralArithmetic T = double, std::floating_point U = T>
using Function1Proxy = pro::proxy<detail::Function1Facader<T, U>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeFunction1Proxy(T function1)
    -> Function1Proxy<typename T::value_type, typename T::param_type> {
    return pro::make_proxy<
        detail::Function1Facader<typename T::value_type, typename T::param_type>>(
        std::move(function1)
    );
}

} // namespace boyle::math
