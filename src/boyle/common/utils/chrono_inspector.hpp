/**
 * @file chrono_inspector.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-03
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <chrono>
#include <source_location>
#include <string_view>

#include "spdlog/fmt/chrono.h"

#include "boyle/common/utils/logging.hpp"

namespace boyle::common {

template <typename Clock = std::chrono::steady_clock, typename Duration = Clock::duration>
    requires std::chrono::is_clock_v<Clock>
class [[nodiscard]] ChronoInspector final {
  public:
    using clock_type = Clock;
    using time_point = typename clock_type::time_point;
    using duration = Duration;

    ChronoInspector() noexcept = delete;
    ChronoInspector(const ChronoInspector& other) noexcept = delete;
    auto operator=(const ChronoInspector& other) noexcept -> ChronoInspector& = delete;
    ChronoInspector(ChronoInspector&& other) noexcept = delete;
    auto operator=(ChronoInspector&& other) noexcept -> ChronoInspector& = delete;
    ~ChronoInspector() noexcept {
        m_logger->log(
            m_source_loc, boyle::common::LogLevel::trace, "{0:s}: {1}.", m_info, elapsed()
        );
    }

    [[using gnu: always_inline]]
    explicit ChronoInspector(
        std::string_view info,
        const std::source_location& source_loc = std::source_location::current()
    ) noexcept
        : m_logger{::boyle::common::getDefaultLogger()},
          m_source_loc{
              source_loc.file_name(), static_cast<int>(source_loc.line()),
              source_loc.function_name()
          },
          m_info{info}, m_start{clock_type::now()} {}

    [[using gnu: always_inline]]
    explicit ChronoInspector(
        std::shared_ptr<boyle::common::Logger> logger, std::string_view info,
        const std::source_location& source_loc = std::source_location::current()
    ) noexcept
        : m_logger{std::move(logger)},
          m_source_loc{
              source_loc.file_name(), static_cast<int>(source_loc.line()),
              source_loc.function_name()
          },
          m_info{info}, m_start{clock_type::now()} {}

    [[using gnu: always_inline]]
    auto reset() noexcept -> void {
        m_start = clock_type::now();
        return;
    }

    template <typename OtherDuration = duration>
    [[using gnu: pure, always_inline]]
    auto elapsed() const noexcept -> OtherDuration {
        return std::chrono::duration_cast<OtherDuration>(clock_type::now() - m_start);
    }

  private:
    std::shared_ptr<boyle::common::Logger> m_logger;
    spdlog::source_loc m_source_loc;
    std::string_view m_info;
    time_point m_start;
};

} // namespace boyle::common
