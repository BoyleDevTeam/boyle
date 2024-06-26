/**
 * @file type_traits.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-12-28
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <type_traits>

namespace boyle::math {

namespace detail {

template <template <typename...> class Base, typename Derived>
struct isBaseOfTemplateImpl final {
    template <typename... Ts>
    static consteval auto test(const Base<Ts...>*) noexcept -> std::true_type;
    static consteval auto test(...) noexcept -> std::false_type;
    using type = decltype(test(std::declval<Derived*>()));
};

} // namespace detail

template <template <typename...> class Base, typename Derived>
using isBaseOfTemplate = typename detail::isBaseOfTemplateImpl<Base, Derived>::type;

template <template <typename...> class Base, typename Derived>
inline constexpr bool isBaseOfTemplateV = isBaseOfTemplate<Base, Derived>::value;

template <typename Instance, template <typename...> class TemplateType>
struct isInstanceOfTemplate final : public std::false_type {};

template <typename Instance, template <typename...> class TemplateType>
inline constexpr bool isInstanceOfTemplateV = isInstanceOfTemplate<Instance, TemplateType>::value;

template <template <typename...> class TemplateType, typename... Argv>
struct isInstanceOfTemplate<TemplateType<Argv...>, TemplateType> : public std::true_type {};

template <typename T, typename = void>
struct isIterable final : public std::false_type {};

template <typename T>
inline constexpr bool IsIterableV = isIterable<T>::value;

template <typename T>
struct isIterable<
    T, std::void_t<decltype(std::declval<T>().begin()), decltype(std::declval<T>().end())>>
    final : public std::true_type {};

template <typename T, typename = void>
struct isVecArithmetic final : public std::false_type {};

template <typename T>
inline constexpr bool isVecArithmeticV = isVecArithmetic<T>::value;

template <typename T>
struct isVecArithmetic<
    T, std::void_t<
           decltype(std::declval<T>() + std::declval<T>()),
           decltype(std::declval<T>() - std::declval<T>()),
           decltype(std::declval<T>().dot(std::declval<T>())),
           decltype(std::declval<T>() += std::declval<T>()),
           decltype(std::declval<T>() -= std::declval<T>()),
           decltype(std::declval<T>().euclideanTo(std::declval<T>())),
           decltype(std::declval<T>() * std::declval<double>()),
           decltype(std::declval<double>() * std::declval<T>()),
           decltype(std::declval<T>() / std::declval<double>()),
           decltype(std::declval<T>() *= std::declval<double>()),
           decltype(std::declval<T>() /= std::declval<double>())>>
    final : public std::true_type {};

template <typename T>
struct isGeneralVecArithmetic final
    : public std::integral_constant<
          bool, std::is_arithmetic<T>::value || isVecArithmetic<T>::value> {};

template <typename T>
inline constexpr bool isGeneralVecArithmeticV = isGeneralVecArithmetic<T>::value;

} // namespace boyle::math
