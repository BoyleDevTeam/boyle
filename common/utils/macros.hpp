/**
 * @file macros.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-19
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#define ENABLE_COPY(ClassName)                            \
    ClassName(const ClassName& other) noexcept = default; \
    ClassName& operator=(const ClassName& other) noexcept = default

#define DISABLE_COPY(ClassName)                          \
    ClassName(const ClassName& other) noexcept = delete; \
    ClassName& operator=(const ClassName& other) noexcept = delete

#define ENABLE_MOVE(ClassName)                       \
    ClassName(ClassName&& other) noexcept = default; \
    ClassName& operator=(ClassName&& other) noexcept = default

#define DISABLE_MOVE(ClassName)                     \
    ClassName(ClassName&& other) noexcept = delete; \
    ClassName& operator=(ClassName&& other) noexcept = delete

#define ENABLE_COPY_AND_MOVE(ClassName) \
    ENABLE_COPY(ClassName);             \
    ENABLE_MOVE(ClassName)

#define DISABLE_COPY_AND_MOVE(ClassName) \
    DISABLE_COPY(ClassName);             \
    DISABLE_MOVE(ClassName)

#define ENABLE_IMPLICIT_CONSTRUCTORS(ClassName) \
    ClassName() noexcept = default;             \
    ENABLE_COPY_AND_MOVE(ClassName)

#define DISABLE_IMPLICIT_CONSTRUCTORS(ClassName) \
    ClassName() noexcept = delete;               \
    DISABLE_COPY_AND_MOVE(ClassName)

#define MAKE_SINGLETON(ClassName)              \
    DISABLE_COPY_AND_MOVE(ClassName);          \
    static ClassName* getInstance() noexcept { \
        static ClassName instance;             \
        return &instance;                      \
    }
