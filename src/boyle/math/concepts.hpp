/**
 * @file concepts.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-11-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <type_traits>

#include "boyle/math/type_traits.hpp"

namespace boyle::math {

template <typename Derived, template <typename...> typename Base>
concept derivedFromTemplate = isBaseOfTemplateV<Base, Derived>;

template <typename T>
concept Iterable = requires(T a) {
    { a.cbegin() } -> std::same_as<typename T::const_iterator>;
    { a.cend() } -> std::same_as<typename T::const_iterator>;
};

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept VecArithmetic = (
    std::floating_point<typename T::value_type> &&

    requires(T a, T b) {
        { a + b } -> std::same_as<T>;
        { a - b } -> std::same_as<T>;
        { a.dot(b) } -> std::same_as<typename T::value_type>;
        { a += b } -> std::same_as<T&>;
        { a -= b } -> std::same_as<T&>;
        { a.distanceTo(b) } -> std::same_as<typename T::value_type>;
    } &&

    requires(T a, typename T::value_type scale) {
        { a* scale } -> std::same_as<T>;
        { scale* a } -> std::same_as<T>;
        { a / scale } -> std::same_as<T>;
        { a *= scale } -> std::same_as<T&>;
        { a /= scale } -> std::same_as<T&>;
    }
);

template <typename T>
concept GeneralArithmetic = Arithmetic<T> || VecArithmetic<T>;

} // namespace boyle::math
