/**
 * @file dualism.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-21
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

namespace tiny_pnc {
namespace common {

enum class Stiffness : unsigned int {
    HARD,
    SOFT
};

enum class Chirality : unsigned int {
    LEFT,
    RIGHT
};

enum class Actio : unsigned int {
    BLOCKING,
    PUSHING
};

} // namespace common
} // namespace tiny_pnc
