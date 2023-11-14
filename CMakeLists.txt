cmake_minimum_required(VERSION 3.26.4)

project(tiny_pnc
  LANGUAGES
    C
    CXX
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_DEBUG "-Og")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-Og")

set(CMAKE_COMPILE_WARNING_AS_ERROR True)

add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Wno-uninitialized
    -Wno-unused-parameter
)

set(CPM_SOURCE_CACHE third_party)
set(CPM_USE_LOCAL_PACKAGES True)

add_compile_definitions(
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE
    DOCTEST_CONFIG_SUPER_FAST_ASSERTS
)

include(CTest)
enable_testing()

include(cmake/CPM.cmake)
include(cmake/third_party.cmake)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(common)
add_subdirectory(math)
