/**
 * @file road_graph_data.h
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
#include <type_traits>
#include <unordered_map>

#include "common/road_elements/road_element.hpp"
#include "common/utils/macros.hpp"

namespace tiny_pnc {
namespace common {

class RoadGraphData final {
  public:
    explicit RoadGraphData(std::unordered_map<std::uint64_t, RoadElement*> road_elements) noexcept
        : road_elements_(std::move(road_elements)) {}
    RoadGraphData() noexcept = default;
    DISABLE_COPY_AND_MOVE(RoadGraphData);
    ~RoadGraphData() noexcept {
        for (auto& [id, road_element] : road_elements_) {
            delete road_element;
            road_element = nullptr;
        }
    }
    RoadElement* getRoadElementById(std::uint64_t id) { return road_elements_.at(id); }
    const RoadElement* getRoadElementById(std::uint64_t id) const { return road_elements_.at(id); }

  private:
    std::unordered_map<std::uint64_t, RoadElement*> road_elements_;
};

} // namespace common
} // namespace tiny_pnc
