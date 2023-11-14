/**
 * @file tags.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-12
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

namespace tiny_pnc {
namespace common {

struct allocate_tag final {
    constexpr explicit allocate_tag() noexcept = default;
};

struct raw_construct_tag final {
    constexpr explicit raw_construct_tag() noexcept = default;
};

} // namespace common
} // namespace tiny_pnc
