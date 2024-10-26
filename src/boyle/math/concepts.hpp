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

template <class Derived, template <typename...> class Base>
concept derivedFromTemplate = isBaseOfTemplateV<Base, Derived>;

template <class Instance, template <typename...> class TemplateType>
concept InstanceOfTemplate = isInstanceOfTemplateV<Instance, TemplateType>;

template <typename T>
concept Iterable = requires(T a) {
    { a.begin() } -> std::same_as<typename T::iterator>;
    { a.end() } -> std::same_as<typename T::iterator>;
};

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept ComplexArithmetic = (
    std::floating_point<typename T::value_type> &&

    requires(T a) {
        {a.real()} -> std::same_as<typename T::value_type>;
        {a.imag()} -> std::same_as<typename T::value_type>;
    } &&

    requires(T a, T b) {
        {a + b} -> std::same_as<T>;
        {a - b} -> std::same_as<T>;
        {a * b} -> std::same_as<T>;
        {a / b} -> std::same_as<T>;
        {a += b} -> std::same_as<T&>;
        {a -= b} -> std::same_as<T&>;
        {a *= b} -> std::same_as<T&>;
        {a /= b} -> std::same_as<T&>;
    } &&

    requires(T a, typename T::value_type b) {
        {a + b} -> std::same_as<T>;
        {a - b} -> std::same_as<T>;
        {a * b} -> std::same_as<T>;
        {a / b} -> std::same_as<T>;
        {a += b} -> std::same_as<T&>;
        {a -= b} -> std::same_as<T&>;
        {a *= b} -> std::same_as<T&>;
        {a /= b} -> std::same_as<T&>;
        {b + a} -> std::same_as<T>;
        {b - a} -> std::same_as<T>;
        {b * a} -> std::same_as<T>;
        {b / a} -> std::same_as<T>;
    }
);

template <typename T>
concept VecArithmetic = (
    // clang-format off-next-line
    (std::floating_point<typename T::value_type> || ComplexArithmetic<typename T::value_type>) &&

    requires(T a, T b) {
        { a + b } -> std::same_as<T>;
        { a - b } -> std::same_as<T>;
        { a.dot(b) } -> std::same_as<typename T::value_type>;
        { a += b } -> std::same_as<T&>;
        { a -= b } -> std::same_as<T&>;
        { a.euclideanTo(b) } -> std::same_as<typename T::value_type>;
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
concept GeneralArithmetic = std::floating_point<T> || ComplexArithmetic<T> || VecArithmetic<T>;

} // namespace boyle::math
