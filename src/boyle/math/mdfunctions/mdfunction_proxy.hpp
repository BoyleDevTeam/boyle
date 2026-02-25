/**
 * @file mdfunction_proxy.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-09-07
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <utility>

#include "proxy/proxy.h"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_num_dimensions, num_dimensions);
PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_gradient, gradient);
PRO_DEF_MEM_DISPATCH(c_has_extrema, has_extrema);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <VecArithmetic T>
class MdFunctionFacade final : public pro::facade_builder
    ::support_copy<pro::constraint_level::nontrivial>
    ::add_convention<c_num_dimensions, auto() const noexcept -> typename T::size_type>
    ::template add_convention<pro::operator_dispatch<"()">, auto(const T&) const noexcept -> typename T::value_type>
    ::template add_convention<c_eval, auto(const T&) const noexcept -> typename T::value_type>
    ::template add_convention<c_gradient, auto(const T&) const noexcept -> T>
    ::template add_convention<c_gradient, auto(const T&, std::size_t) const noexcept -> typename T::value_type>
    ::template add_convention<c_has_extrema, auto(const T&) const noexcept -> bool>
    ::build {};

// clang-format on

} // namespace detail

template <VecArithmetic T>
using MdFunctionProxy = pro::proxy<detail::MdFunctionFacade<T>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeMdFunctionProxy(T&& mdfunction)
    -> MdFunctionProxy<typename std::remove_cvref_t<T>::param_type> {
    return pro::make_proxy<detail::MdFunctionFacade<typename std::remove_cvref_t<T>::param_type>>(
        std::forward<T>(mdfunction)
    );
}

} // namespace boyle::math
