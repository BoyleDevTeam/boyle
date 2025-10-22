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

template <typename T>
concept Iterable = requires(T a) {
    { a.begin() } -> std::same_as<typename std::remove_cvref_t<T>::iterator>;
    { a.end() } -> std::same_as<typename std::remove_cvref_t<T>::iterator>;
};

template <typename T>
concept Arithmetic = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template <typename T>
concept ComplexArithmetic = (
    requires(std::remove_cvref_t<T> a) {
        { a.real() } -> std::same_as<typename std::remove_cvref_t<T>::value_type>;
        { a.imag() } -> std::same_as<typename std::remove_cvref_t<T>::value_type>;
    } &&

    requires(std::remove_cvref_t<T> a, std::remove_cvref_t<T> b) {
        { a + b } -> std::same_as<std::remove_cvref_t<T>>;
        { a - b } -> std::same_as<std::remove_cvref_t<T>>;
        { a * b } -> std::same_as<std::remove_cvref_t<T>>;
        { a / b } -> std::same_as<std::remove_cvref_t<T>>;
    } &&

    requires(std::remove_cvref_t<T> a, typename std::remove_cvref_t<T>::value_type b) {
        { a + b } -> std::same_as<std::remove_cvref_t<T>>;
        { a - b } -> std::same_as<std::remove_cvref_t<T>>;
        { a * b } -> std::same_as<std::remove_cvref_t<T>>;
        { a / b } -> std::same_as<std::remove_cvref_t<T>>;
        { b + a } -> std::same_as<std::remove_cvref_t<T>>;
        { b - a } -> std::same_as<std::remove_cvref_t<T>>;
        { b * a } -> std::same_as<std::remove_cvref_t<T>>;
        { b / a } -> std::same_as<std::remove_cvref_t<T>>;
    }
);

template <typename T>
concept VecArithmetic = (
    requires(std::remove_cvref_t<T> a) {
        { a.data() };
        { a.size() } -> std::same_as<typename std::remove_cvref_t<T>::size_type>;
        { a.stride() } -> std::same_as<typename std::remove_cvref_t<T>::size_type>;
    } &&

    requires(std::remove_cvref_t<T> a, typename std::remove_cvref_t<T>::size_type i) {
        { a[i] } -> std::same_as<typename std::remove_cvref_t<T>::value_type&>;
        { a.coeff(i) } -> std::same_as<typename std::remove_cvref_t<T>::value_type>;
    } &&

    requires(std::remove_cvref_t<T> a, std::remove_cvref_t<T> b) {
        { a + b } -> std::same_as<std::remove_cvref_t<T>>;
        { a - b } -> std::same_as<std::remove_cvref_t<T>>;
        { a.dot(b) } -> std::same_as<typename std::remove_cvref_t<T>::value_type>;
        { a.euclideanTo(b) };
    } &&

    requires(std::remove_cvref_t<T> a, typename std::remove_cvref_t<T>::value_type b) {
        { a * b } -> std::same_as<std::remove_cvref_t<T>>;
        { b * a } -> std::same_as<std::remove_cvref_t<T>>;
        { a / b } -> std::same_as<std::remove_cvref_t<T>>;
    }
);

template <typename T>
concept MatArithmetic = (
    requires(std::remove_cvref_t<T> a) {
        { a.data() };
        { a.nrows() } -> std::same_as<typename std::remove_cvref_t<T>::size_type>;
        { a.ncols() } -> std::same_as<typename std::remove_cvref_t<T>::size_type>;
        { a.stride() } -> std::same_as<typename std::remove_cvref_t<T>::size_type>;
    } &&

    requires(std::remove_cvref_t<T> a, typename std::remove_cvref_t<T>::size_type i, typename std::remove_cvref_t<T>::size_type j) {
        { a[i, j] } -> std::same_as<typename std::remove_cvref_t<T>::value_type&>;
        { a.coeff(i, j) } -> std::same_as<typename std::remove_cvref_t<T>::value_type>;
    } &&

    requires(std::remove_cvref_t<T> a, std::remove_cvref_t<T> b) {
        { a + b } -> std::same_as<std::remove_cvref_t<T>>;
        { a - b } -> std::same_as<std::remove_cvref_t<T>>;
        { a.dot(b) };
    } &&

    requires(std::remove_cvref_t<T> a, typename std::remove_cvref_t<T>::value_type b) {
        { a * b } -> std::same_as<std::remove_cvref_t<T>>;
        { b * a } -> std::same_as<std::remove_cvref_t<T>>;
        { a / b } -> std::same_as<std::remove_cvref_t<T>>;
    }
);

template <typename T>
concept ScalarArithmetic = Arithmetic<T> || ComplexArithmetic<T>;

template <typename T>
concept GeneralArithmetic = ScalarArithmetic<T> || VecArithmetic<T> || MatArithmetic<T>;

} // namespace boyle::math
