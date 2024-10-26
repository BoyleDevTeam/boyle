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

PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_derivative, derivative);
PRO_DEF_MEM_DISPATCH(c_integral, integral);
PRO_DEF_MEM_DISPATCH(c_minT, minT);
PRO_DEF_MEM_DISPATCH(c_maxT, maxT);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <::boyle::math::GeneralArithmetic T, std::floating_point U = T>
class [[nodiscard]] Function1Facade final : public pro::facade_builder
    ::template add_convention<c_eval, auto(U) const noexcept -> T>
    ::template add_convention<c_derivative, auto(U) const noexcept -> T>
    ::template add_convention<c_integral, auto(U, U) const noexcept -> T>
    ::template add_convention<c_minT, auto() const noexcept -> U>
    ::template add_convention<c_maxT, auto() const noexcept -> U>
    ::build {};

// clang-format on

} // namespace detail

template <GeneralArithmetic T = double, std::floating_point U = T>
using Function1Proxy = pro::proxy<detail::Function1Facade<T, U>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeFunction1Proxy(T function1
) -> Function1Proxy<typename T::value_type, typename T::param_type> {
    return pro::make_proxy<detail::Function1Facade<typename T::value_type, typename T::param_type>>(
        std::move(function1)
    );
}

} // namespace boyle::math
