/**
 * @file logging.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-05-09
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include "spdlog/async.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#define BOYLE_LOG_LEVEL_TRACE SPDLOG_LEVEL_TRACE
#define BOYLE_LOG_LEVEL_DEBUG SPDLOG_LEVEL_DEBUG
#define BOYLE_LOG_LEVEL_INFO SPDLOG_LEVEL_INFO
#define BOYLE_LOG_LEVEL_WARN SPDLOG_LEVEL_WARN
#define BOYLE_LOG_LEVEL_ERROR SPDLOG_LEVEL_ERROR
#define BOYLE_LOG_LEVEL_CRITICAL SPDLOG_LEVEL_CRITICAL
#define BOYLE_LOG_LEVEL_OFF SPDLOG_LEVEL_OFF

#if !defined(BOYLE_LOG_ACTIVE_LEVEL)
#define BOYLE_LOG_ACTIVE_LEVEL BOYLE_LOG_LEVEL_INFO
#endif

#ifndef BOYLE_LOG_NO_SOURCE_LOC
#define BOYLE_LOG_LOGGER_CALL(logger, level, ...) \
    (logger)->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, level, __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_CALL(logger, level, ...) \
    (logger)->log(spdlog::source_loc{}, level, __VA_ARGS__)
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_TRACE
#define BOYLE_LOG_LOGGER_TRACE(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::trace, __VA_ARGS__)
#define BOYLE_LOG_TRACE(...) BOYLE_LOG_LOGGER_TRACE(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_TRACE(logger, ...) (void)0
#define BOYLE_LOG_TRACE(...) (void)0
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_DEBUG
#define BOYLE_LOG_LOGGER_DEBUG(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::debug, __VA_ARGS__)
#define BOYLE_LOG_DEBUG(...) BOYLE_LOG_LOGGER_DEBUG(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_DEBUG(logger, ...) (void)0
#define BOYLE_LOG_DEBUG(...) (void)0
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_INFO
#define BOYLE_LOG_LOGGER_INFO(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::info, __VA_ARGS__)
#define BOYLE_LOG_INFO(...) BOYLE_LOG_LOGGER_INFO(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_INFO(logger, ...) (void)0
#define BOYLE_LOG_INFO(...) (void)0
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_WARN
#define BOYLE_LOG_LOGGER_WARN(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::warn, __VA_ARGS__)
#define BOYLE_LOG_WARN(...) BOYLE_LOG_LOGGER_WARN(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_WARN(logger, ...) (void)0
#define BOYLE_LOG_WARN(...) (void)0
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_ERROR
#define BOYLE_LOG_LOGGER_ERROR(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::err, __VA_ARGS__)
#define BOYLE_LOG_ERROR(...) BOYLE_LOG_LOGGER_ERROR(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_ERROR(logger, ...) (void)0
#define BOYLE_LOG_ERROR(...) (void)0
#endif

#if BOYLE_LOG_ACTIVE_LEVEL <= BOYLE_LOG_LEVEL_CRITICAL
#define BOYLE_LOG_LOGGER_CRITICAL(logger, ...) \
    BOYLE_LOG_LOGGER_CALL(logger, spdlog::level::critical, __VA_ARGS__)
#define BOYLE_LOG_CRITICAL(...) BOYLE_LOG_LOGGER_CRITICAL(spdlog::default_logger_raw(), __VA_ARGS__)
#else
#define BOYLE_LOG_LOGGER_CRITICAL(logger, ...) (void)0
#define BOYLE_LOG_CRITICAL(...) (void)0
#endif

#define BOYLE_LOG_LOGGER_CALL_IF(logger, level, condition, ...) \
    (condition) ? BOYLE_LOG_LOGGER_CALL(logger, level, __VA_ARGS__) : (void)0

#define BOYLE_LOG_LOGGER_TRACE_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::trace, condition, __VA_ARGS__)

#define BOYLE_LOG_LOGGER_DEBUG_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::debug, condition, __VA_ARGS__)

#define BOYLE_LOG_LOGGER_INFO_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::info, condition, __VA_ARGS__)

#define BOYLE_LOG_LOGGER_WARN_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::warn, condition, __VA_ARGS__)

#define BOYLE_LOG_LOGGER_ERROR_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::err, condition, __VA_ARGS__)

#define BOYLE_LOG_LOGGER_CRITICAL_IF(logger, condition, ...) \
    BOYLE_LOG_LOGGER_CALL_IF(logger, spdlog::level::critical, condition, __VA_ARGS__)

#define BOYLE_LOG_TRACE_IF(condition, ...) \
    BOYLE_LOG_LOGGER_TRACE_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

#define BOYLE_LOG_DEBUG_IF(condition, ...) \
    BOYLE_LOG_LOGGER_DEBUG_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

#define BOYLE_LOG_INFO_IF(condition, ...) \
    BOYLE_LOG_LOGGER_INFO_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

#define BOYLE_LOG_WARN_IF(condition, ...) \
    BOYLE_LOG_LOGGER_WARN_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

#define BOYLE_LOG_ERROR_IF(condition, ...) \
    BOYLE_LOG_LOGGER_ERROR_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

#define BOYLE_LOG_CRITICAL_IF(condition, ...) \
    BOYLE_LOG_LOGGER_CRITICAL_IF(spdlog::default_logger_raw(), condition, __VA_ARGS__)

namespace boyle::common {

using Logger = spdlog::logger;
using LogLevel = spdlog::level::level_enum;

inline auto initLogger(const std::string& name, const std::string& log_file = "") -> void {
    std::vector<std::shared_ptr<spdlog::sinks::sink>> sinks{};

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::trace);
    console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] [thread %t] %@ %v");
    sinks.push_back(console_sink);

    if (!log_file.empty()) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
        file_sink->set_level(spdlog::level::trace);
        file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] [thread %t] %@ %v");
        sinks.push_back(file_sink);
    }

    spdlog::init_thread_pool(8192, 1);
    auto logger = std::make_shared<spdlog::async_logger>(
        name, sinks.cbegin(), sinks.cend(), spdlog::thread_pool()
    );
    spdlog::drop_all();
    spdlog::register_logger(logger);
    spdlog::set_default_logger(logger);

    return;
}

inline auto getLogger(const std::string& name) -> std::shared_ptr<Logger> {
    return spdlog::get(name);
}

} // namespace boyle::common
