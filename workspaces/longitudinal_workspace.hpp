/**
 * @file longitudinal_workspace.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-02
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
#include "workspaces/workspace.hpp"

namespace tiny_pnc {

class LongitudinalWorkspace final : public Workspace {
  public:
    LongitudinalWorkspace() noexcept = default;
    explicit LongitudinalWorkspace(std::unordered_map<std::string, Workspace*> sub_workspace_map
    ) noexcept
        : Workspace(std::move(sub_workspace_map)) {}
    DISABLE_COPY_AND_MOVE(LongitudinalWorkspace);
    ~LongitudinalWorkspace() noexcept override = default;

  private:
    /* Put function object here */
};

} // namespace tiny_pnc
