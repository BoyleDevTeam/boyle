/**
 * @file dualism.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-21
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

namespace boyle::kinetics {

enum class Chirality : unsigned int {
    LEFT,
    RIGHT
};

enum class Actio : unsigned int {
    BLOCKING,
    PUSHING
};

} // namespace boyle::kinetics
