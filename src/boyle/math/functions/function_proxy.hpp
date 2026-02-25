/**
 * @file function_proxy.hpp
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

#include <utility>

#include "proxy/proxy.h"

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/detail/dense_degenerate_trait.hpp"

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

template <::boyle::math::GeneralArithmetic T>
class FunctionFacade final : public pro::facade_builder
    ::support_copy<pro::constraint_level::nontrivial>
    ::template add_convention<pro::operator_dispatch<"()">, auto(::boyle::math::detail::DenseDegenerateTraitT<T>) const noexcept -> T>
    ::template add_convention<c_eval, auto(::boyle::math::detail::DenseDegenerateTraitT<T>) const noexcept -> T>
    ::template add_convention<c_derivative, auto(::boyle::math::detail::DenseDegenerateTraitT<T>) const noexcept -> T>
    ::template add_convention<c_integral, auto(::boyle::math::detail::DenseDegenerateTraitT<T>, ::boyle::math::detail::DenseDegenerateTraitT<T>) const noexcept -> T>
    ::template add_convention<c_minT, auto() const noexcept -> ::boyle::math::detail::DenseDegenerateTraitT<T>>
    ::template add_convention<c_maxT, auto() const noexcept -> ::boyle::math::detail::DenseDegenerateTraitT<T>>
    ::build {};

// clang-format on

} // namespace detail

template <GeneralArithmetic T>
using FunctionProxy = pro::proxy<detail::FunctionFacade<T>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeFunctionProxy(T&& function1)
    -> FunctionProxy<typename std::remove_cvref_t<T>::value_type> {
    return pro::make_proxy<detail::FunctionFacade<typename std::remove_cvref_t<T>::value_type>>(
        std::forward<T>(function1)
    );
}

} // namespace boyle::math
