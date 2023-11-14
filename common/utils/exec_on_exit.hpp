/**
 * @file exec_on_exit.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <functional>

#include "common/utils/macros.hpp"

namespace tiny_pnc {
namespace common {

class [[nodiscard]] ExecOnExit final {
  public:
    using Callback = std::function<void()>;
    ExecOnExit(Callback callback) noexcept : callback_(callback) {}
    DISABLE_IMPLICIT_CONSTRUCTORS(ExecOnExit);
    ~ExecOnExit() { callback_(); }

  private:
    Callback callback_;
};

} // namespace common
} // namespace tiny_pnc
