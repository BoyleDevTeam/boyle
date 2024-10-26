/**
 * @file mdfunction_proxy.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-14
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "proxy.h"

namespace boyle::math {

namespace detail {

// NOLINTBEGIN(modernize-use-trailing-return-type)

PRO_DEF_MEM_DISPATCH(c_num_dimensions, n_num_dimensions);
PRO_DEF_MEM_DISPATCH(c_eval, eval);
PRO_DEF_MEM_DISPATCH(c_gradient, gradient);
PRO_DEF_MEM_DISPATCH(c_has_extrema, has_extrema);

// NOLINTEND(modernize-use-trailing-return-type)

// clang-format off

template <std::floating_point T = double, typename U = std::span<const T>>
class [[nodiscard]] MdFunctionFacade final : public pro::facade_builder
    ::add_convention<c_num_dimensions, auto() const noexcept -> std::size_t>
    ::template add_convention<c_eval, auto(U) const noexcept -> T>
    ::template add_convention<c_gradient, auto(U) const noexcept -> std::vector<T>>
    ::template add_convention<c_gradient, auto(U, std::size_t) const noexcept -> T>
    ::template add_convention<c_has_extrema, auto(U) const noexcept -> bool>
    ::build {};

// clang-format on

} // namespace detail

template <std::floating_point T = double, typename U = std::span<const T>>
using MdFunctionProxy = pro::proxy<detail::MdFunctionFacade<T, U>>;

template <typename T>
[[using gnu: always_inline]]
inline auto makeMdFunctionProxy(T&& mdfunction)
    -> MdFunctionProxy<
        typename std::remove_reference_t<T>::value_type,
        typename std::remove_reference_t<T>::param_type> {
    return pro::make_proxy<detail::MdFunctionFacade<
        typename std::remove_reference_t<T>::value_type,
        typename std::remove_reference_t<T>::param_type>>(std::forward(mdfunction));
}

} // namespace boyle::math
