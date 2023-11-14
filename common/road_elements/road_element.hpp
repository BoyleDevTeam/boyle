/**
 * @file road_element.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-20
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <cstdint>

#include "common/utils/macros.hpp"

namespace tiny_pnc {
namespace common {

class RoadElement {
  public:
    explicit RoadElement(std::uint64_t id) noexcept : id_(id) {}
    ENABLE_IMPLICIT_CONSTRUCTORS(RoadElement);
    virtual ~RoadElement() noexcept = default;
    const std::uint64_t& id() const noexcept { return id_; }

  protected:
    std::uint64_t id_;
};

} // namespace common
} // namespace tiny_pnc
