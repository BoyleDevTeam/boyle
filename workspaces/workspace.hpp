/**
 * @file workspace.hpp
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

namespace tiny_pnc {

class Workspace {
  public:
    Workspace() noexcept = default;
    explicit Workspace(std::unordered_map<std::string, Workspace*> sub_workspace_map) noexcept
        : sub_workspace_map_(std::move(sub_workspace_map)) {}
    DISABLE_COPY_AND_MOVE(Workspace);
    virtual ~Workspace() noexcept {
        for (auto& [key, sub_workspace] : sub_workspace_map_) {
            delete sub_workspace;
            sub_workspace = nullptr;
        }
    }
    Workspace* getSubWorkspace(const std::string& key) { return sub_workspace_map_.at(key); }
    const Workspace* getSubWorkspace(const std::string& key) const {
        return sub_workspace_map_.at(key);
    }

  private:
    std::unordered_map<std::string, Workspace*> sub_workspace_map_;
};

} // namespace tiny_pnc
