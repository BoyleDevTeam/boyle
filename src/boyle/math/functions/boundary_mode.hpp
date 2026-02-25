/**
 * @file boundary_mode.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-11-20
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include "boyle/math/concepts.hpp"

namespace boyle::math {

template <GeneralArithmetic T>
struct BoundaryMode final {
    unsigned int order;
    T derivative;
};

} // namespace boyle::math
