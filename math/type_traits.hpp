/**
 * @file type_traits.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-11
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

namespace tiny_pnc {
namespace math {

template <typename T>
struct isVecArithmetic final {
    static constexpr bool value = false;
};

template <typename T>
inline constexpr bool isVecArithmeticV = isVecArithmetic<T>::value;

template <typename T>
struct isFunction final {
    static constexpr bool value = false;
};

template <typename T>
inline constexpr bool isFunctionV = isFunction<T>::value;

template <typename T>
struct isCurve final {
    static constexpr bool value = false;
};

template <typename T>
inline constexpr bool isCurveV = isCurve<T>::value;

} // namespace math
} // namespace tiny_pnc
