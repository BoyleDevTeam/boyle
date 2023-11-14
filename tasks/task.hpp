/**
 * @file task.hpp
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
#include <utility>

#include "common/utils/macros.hpp"
#include "workspaces/workspace.hpp"

namespace tiny_pnc {

class Task {
  public:
    Task() noexcept = default;
    explicit Task(std::unordered_map<std::string, Task*> sub_task_map) noexcept
        : sub_task_map_(std::move(sub_task_map)) {}
    DISABLE_COPY_AND_MOVE(Task);
    virtual ~Task() noexcept {
        for (auto& [key, subtask] : sub_task_map_) {
            delete subtask;
            subtask = nullptr;
        }
    }
    virtual void prepare(std::initializer_list<const Workspace*> workspace_list) noexcept = 0;
    virtual void process() = 0;
    virtual void populate(std::initializer_list<Workspace*> workspace_list) const noexcept = 0;

  protected:
    Task* getSubTask(const std::string& key) { return sub_task_map_.at(key); }
    const Task* getSubTask(const std::string& key) const { return sub_task_map_.at(key); }

  private:
    std::unordered_map<std::string, Task*> sub_task_map_;
};

} // namespace tiny_pnc
