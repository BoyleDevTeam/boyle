/**
 * @file exec_on_exit.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <functional>
#include <utility>

#include "boyle/common/utils/macros.hpp"

namespace boyle::common {

class [[nodiscard]] ExecOnExit final {
  public:
    using Callback = std::function<auto()->void>;
    ExecOnExit(Callback callback) noexcept : m_callback{std::move(callback)} {}
    DISABLE_IMPLICIT_CONSTRUCTORS(ExecOnExit);
    ~ExecOnExit() { m_callback(); }

  private:
    Callback m_callback;
};

} // namespace boyle::common
