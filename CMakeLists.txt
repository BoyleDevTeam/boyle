cmake_minimum_required(VERSION 3.31.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS True CACHE INTERNAL "Export compile_commands.json")

project(Boyle
  LANGUAGES
    C
    CXX
  VERSION
    0.1.0
)
set(BOYLE_SOVERSION 0)

set(CMAKE_C_STANDARD "23" CACHE INTERNAL "C standard")
set(CMAKE_CXX_STANDARD "23" CACHE INTERNAL "C++ standard")
set(CMAKE_COMPILE_WARNING_AS_ERROR True CACHE INTERNAL "Treat warnings as errors")

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/cmake/toolchains/linux-gcc-x64.cmake" CACHE FILEPATH "Toolchain file")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")

option(CMAKE_CXX_EXTENSIONS "Enable C++ extensions" ON)
option(CMAKE_UNITY_BUILD "Enable unity build" OFF)
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(BOYLE_CHECK_PARAMS "Enable parameters checking" OFF)
option(BOYLE_BUILD_TESTING "Enable testing" ON)
option(BOYLE_ENABLE_INSTALL "Enable install" ON)

if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(BOYLE_ENABLE_INSTALL OFF)
endif()

if(BOYLE_CHECK_PARAMS)
  add_compile_definitions(BOYLE_CHECK_PARAMS=1)
else()
  add_compile_definitions(BOYLE_CHECK_PARAMS=0)
endif()

set(CPM_SOURCE_CACHE "third_party")
set(CPM_USE_LOCAL_PACKAGES True)

include(cmake/CPM.cmake)
include(cmake/third_party.cmake)
include(cmake/utils.cmake)

add_subdirectory(src)
add_subdirectory(tests)

include(cmake/export.cmake)

if(EXISTS "${CMAKE_BINARY_DIR}/compile_commands.json")
  file(CREATE_LINK
    "${CMAKE_BINARY_DIR}/compile_commands.json"
    "${CMAKE_SOURCE_DIR}/compile_commands.json"
    SYMBOLIC
  )
endif()
