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

template <template <typename...> class Base, class Derived>
struct isBaseOfTemplateImpl final {
    template <typename... Ts>
    static consteval auto test(const Base<Ts...>*) noexcept -> std::true_type;
    static consteval auto test(...) noexcept -> std::false_type;
    using type = decltype(test(std::declval<Derived*>()));
};

} // namespace detail

template <template <typename...> class Base, class Derived>
using isBaseOfTemplate = typename detail::isBaseOfTemplateImpl<Base, Derived>::type;

template <template <typename...> class Base, class Derived>
inline constexpr bool isBaseOfTemplateV = isBaseOfTemplate<Base, Derived>::value;

template <typename T, typename = void>
struct isIterable final : public std::false_type {};

template <typename T>
struct isIterable<
    T, std::void_t<decltype(std::declval<T>().begin()), decltype(std::declval<T>().end())>>
    final : public std::true_type {};

template <typename T>
inline constexpr bool IsIterableV = isIterable<T>::value;

template <typename T, typename = void>
struct isComplexArithmetic final : public std::false_type {};

template <typename T>
struct isComplexArithmetic<
    T,
    std::void_t<
        decltype(std::declval<std::remove_cvref_t<T>>().real()),
        decltype(std::declval<std::remove_cvref_t<T>>().imag()),
        decltype(std::declval<std::remove_cvref_t<T>>() + std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() - std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() * std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() / std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() + std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<std::remove_cvref_t<T>>() - std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<std::remove_cvref_t<T>>() * std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<std::remove_cvref_t<T>>() / std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() + std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() - std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() * std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() / std::declval<std::remove_cvref_t<T>>())>>
    final : public std::true_type {};

template <typename T>
inline constexpr bool isComplexArithmeticV = isComplexArithmetic<T>::value;

template <typename T, typename = void>
struct isVecArithmetic final : public std::false_type {};

template <typename T>
struct isVecArithmetic<
    T,
    std::void_t<
        decltype(std::declval<std::remove_cvref_t<T>>().data()),
        decltype(std::declval<std::remove_cvref_t<T>>().size()),
        decltype(std::declval<std::remove_cvref_t<T>>().stride()),
        decltype(std::declval<std::remove_cvref_t<T>>()
                     [std::declval<typename std::remove_cvref_t<T>::size_type>()]),
        decltype(std::declval<std::remove_cvref_t<T>>()
                     .coeff(std::declval<typename std::remove_cvref_t<T>::size_type>())),
        decltype(std::declval<std::remove_cvref_t<T>>() + std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() - std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() * std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<std::remove_cvref_t<T>>() / std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() * std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>()
                     .dot(std::declval<std::remove_cvref_t<T>>())),
        decltype(std::declval<std::remove_cvref_t<T>>()
                     .euclideanTo(std::declval<std::remove_cvref_t<T>>())),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() * std::declval<std::remove_cvref_t<T>>())>>
    final : public std::true_type {};

template <typename T>
inline constexpr bool isVecArithmeticV = isVecArithmetic<T>::value;

template <typename T, typename = void>
struct isMatArithmetic final : public std::false_type {};

template <typename T>
struct isMatArithmetic<
    T,
    std::void_t<
        decltype(std::declval<std::remove_cvref_t<T>>().data()),
        decltype(std::declval<std::remove_cvref_t<T>>().nrows()),
        decltype(std::declval<std::remove_cvref_t<T>>().ncols()),
        decltype(std::declval<std::remove_cvref_t<T>>().stride()),
        decltype(std::declval<std::remove_cvref_t<T>>()
                     [std::declval<typename std::remove_cvref_t<T>::size_type>(),
                      std::declval<typename std::remove_cvref_t<T>::size_type>()]),
        decltype(std::declval<std::remove_cvref_t<T>>().coeff(
            std::declval<typename std::remove_cvref_t<T>::size_type>(),
            std::declval<typename std::remove_cvref_t<T>::size_type>()
        )),
        decltype(std::declval<std::remove_cvref_t<T>>() + std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() - std::declval<std::remove_cvref_t<T>>()),
        decltype(std::declval<std::remove_cvref_t<T>>() * std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<std::remove_cvref_t<T>>() / std::declval<typename std::remove_cvref_t<T>::value_type>()),
        decltype(std::declval<typename std::remove_cvref_t<T>::value_type>() * std::declval<std::remove_cvref_t<T>>())>>
    final : public std::true_type {};

template <typename T>
inline constexpr bool isMatArithmeticV = isMatArithmetic<T>::value;

template <typename T>
struct isScalarArithmetic final
    : public std::integral_constant<
          bool, std::is_arithmetic<T>::value || isComplexArithmetic<T>::value> {};

template <typename T>
inline constexpr bool isScalarArithmeticV = isScalarArithmetic<T>::value;

template <typename T>
struct isGeneralArithmetic final
    : public std::integral_constant<
          bool, isScalarArithmetic<T>::value || isVecArithmetic<T>::value ||
                    isMatArithmetic<T>::value> {};

template <typename T>
inline constexpr bool isGeneralArithmeticV = isGeneralArithmetic<T>::value;

} // namespace boyle::math
