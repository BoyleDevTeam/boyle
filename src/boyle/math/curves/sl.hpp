/**
 * @file sl.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-12-16
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>

namespace boyle::math {

template <std::floating_point T>
struct SlDuplet final {
    using value_type = T;
    value_type s;
    value_type l;
};

using SlDupletf = SlDuplet<float>;
using SlDupletd = SlDuplet<double>;

} // namespace boyle::math
