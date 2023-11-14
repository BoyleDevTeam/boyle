/**
 * @file workspace_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-03
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "workspace.hpp"

#include <memory>
#include <unordered_map>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "workspaces/common_workspace.hpp"

namespace tiny_pnc {

TEST_CASE("CommonWorkspaceAccess") {
    std::unordered_map<std::string, Workspace*> sub_workspace_map;
    sub_workspace_map["common"] = new CommonWorkspace;
    std::shared_ptr<Workspace> root_workspace = std::make_shared<Workspace>(sub_workspace_map);
    CommonWorkspace* rw_common_workspace =
        dynamic_cast<CommonWorkspace*>(root_workspace->getSubWorkspace("common"));
    const CommonWorkspace* ro_common_workspace =
        dynamic_cast<const CommonWorkspace*>(root_workspace->getSubWorkspace("common"));

    CHECK_EQ(rw_common_workspace->vehicle_length, 0.0);
    CHECK_EQ(rw_common_workspace->vehicle_width, 0.0);

    rw_common_workspace->vehicle_length = 6.2;
    rw_common_workspace->vehicle_width = 3.7;

    CHECK_EQ(ro_common_workspace->vehicle_length, 6.2);
    CHECK_EQ(ro_common_workspace->vehicle_width, 3.7);
}

} // namespace tiny_pnc
