cmake_minimum_required(VERSION 3.31.0)

project(Boyle
  VERSION
    0.1.0
  DESCRIPTION
    "The Fundamental Math Basis for Autonomous Driving Vehicles and Robotics."
  HOMEPAGE_URL
    https://github.com/BoyleDevTeam/boyle
  LANGUAGES
    C
    CXX
)
set(BOYLE_SOVERSION 0 CACHE INTERNAL "Boyle shared library version")

set(CMAKE_C_STANDARD 23 CACHE INTERNAL "C standard")
set(CMAKE_CXX_STANDARD 23 CACHE INTERNAL "C++ standard")
set(CMAKE_COMPILE_WARNING_AS_ERROR True CACHE INTERNAL "Treat warnings as errors")

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/cmake/toolchains/linux-gcc-x64.cmake CACHE FILEPATH "Toolchain file")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

option(CMAKE_CXX_EXTENSIONS "Enable C++ extensions" ON)
option(CMAKE_UNITY_BUILD "Enable unity build" OFF)
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(BOYLE_CHECK_PARAMS "Enable parameters checking" OFF)
option(BOYLE_BUILD_TESTING "Enable testing" ON)
option(BOYLE_ENABLE_INSTALL "Enable install" ON)
option(BOYLE_USE_BOOST_UNORDERED "Enable Boost.unordered" ON)

set(BOYLE_USE_SIMD "OFF" CACHE STRING "SIMD implementation to use (AVX512, OFF)")
set_property(CACHE BOYLE_USE_SIMD PROPERTY STRINGS "AVX512;OFF")

set(BOYLE_USE_BLAS_LAPACK "OpenBLAS" CACHE STRING "BLAS/LAPACK library to use (OpenBLAS, Netlib, MKL, OFF)")
set_property(CACHE BOYLE_USE_BLAS_LAPACK PROPERTY STRINGS "OpenBLAS;Netlib;MKL;OFF")

if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(BOYLE_ENABLE_INSTALL OFF)
endif()

if(BOYLE_CHECK_PARAMS)
  add_compile_definitions(BOYLE_CHECK_PARAMS=1)
else()
  add_compile_definitions(BOYLE_CHECK_PARAMS=0)
endif()

set(CPM_SOURCE_CACHE ${CMAKE_SOURCE_DIR}/third_party CACHE PATH "CPM source cache")

include(cmake/CPM.cmake)
include(cmake/third_party.cmake)
include(cmake/utils.cmake)

add_subdirectory(src)
add_subdirectory(tests)

if(BOYLE_ENABLE_INSTALL)
  include(cmake/export.cmake)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS True CACHE INTERNAL "Export compile_commands.json")
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  include(cmake/compile_commands.cmake)
endif()
