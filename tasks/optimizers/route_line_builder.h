/**
 * @file route_line_builder.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-08
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <initializer_list>
#include <string>
#include <unordered_map>

#include "common/utils/macros.hpp"
#include "tasks/task.hpp"
#include "workspaces/workspace.hpp"

namespace tiny_pnc {

class RouteLineBuilder final : public Task {
  public:
    RouteLineBuilder() noexcept = default;
    explicit RouteLineBuilder(std::unordered_map<std::string, Task*> sub_task_map) noexcept;
    DISABLE_COPY_AND_MOVE(RouteLineBuilder);
    ~RouteLineBuilder() noexcept override = default;
    void prepare(std::initializer_list<const Workspace*> workspace_list) noexcept override;
    void process() override;
    void populate(std::initializer_list<Workspace*> workspace_list) const noexcept override;
};

} // namespace tiny_pnc
