/**
 * @file common_workspace.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <string>
#include <unordered_map>
#include <utility>

#include "common/utils/macros.hpp"
#include "math/vec2.hpp"
#include "workspaces/workspace.hpp"

namespace tiny_pnc {

class CommonWorkspace final : public Workspace {
  public:
    CommonWorkspace() noexcept = default;
    explicit CommonWorkspace(std::unordered_map<std::string, Workspace*> sub_workspace_map) noexcept
        : Workspace(std::move(sub_workspace_map)) {}
    explicit CommonWorkspace(
        std::unordered_map<std::string, Workspace*> sub_workspace_map, math::Vec2d c_origin_point,
        double c_vehicle_length, double c_vehicle_width, double c_vehicle_steering_ratio,
        math::Vec2d c_actual_position, math::Vec2d c_actual_heading, double c_actual_velocity,
        double c_actual_acceleration, double c_actual_steering_angle,
        math::Vec2d c_implied_position, math::Vec2d c_implied_heading, double c_implied_velocity,
        double c_implied_acceleration, double c_implied_steering_angle
    ) noexcept
        : Workspace(std::move(sub_workspace_map)), origin_point(c_origin_point),
          vehicle_length(c_vehicle_length), vehicle_width(c_vehicle_width),
          vehicle_steering_ratio(c_vehicle_steering_ratio), actual_position(c_actual_position),
          actual_heading(c_actual_heading), actual_velocity(c_actual_velocity),
          actual_acceleration(c_actual_acceleration),
          actual_steering_angle(c_actual_steering_angle), implied_position(c_implied_position),
          implied_heading(c_implied_heading), implied_velocity(c_implied_velocity),
          implied_acceleration(c_implied_acceleration),
          implied_steering_angle(c_implied_steering_angle) {}
    DISABLE_COPY_AND_MOVE(CommonWorkspace);
    ~CommonWorkspace() noexcept override = default;

    math::Vec2d origin_point;
    double vehicle_length;
    double vehicle_width;
    double vehicle_steering_ratio;
    math::Vec2d actual_position;
    math::Vec2d actual_heading;
    double actual_velocity;
    double actual_acceleration;
    double actual_steering_angle;
    math::Vec2d implied_position;
    math::Vec2d implied_heading;
    double implied_velocity;
    double implied_acceleration;
    double implied_steering_angle;
};

} // namespace tiny_pnc
