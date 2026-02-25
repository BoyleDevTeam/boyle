cmake_minimum_required(VERSION 3.31...4.3)

project(Boyle
  VERSION
    0.1.0
  DESCRIPTION
    "The Fundamental Math Basis for Autonomous Driving Vehicles and Robotics."
  HOMEPAGE_URL
    https://github.com/boyle-org/boyle
  LANGUAGES
    C
    CXX
)
set(BOYLE_SOVERSION 0 CACHE INTERNAL "Boyle shared library version")

set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

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

if(NOT PROJECT_IS_TOP_LEVEL)
  set(BOYLE_ENABLE_INSTALL OFF)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_compile_definitions(BOYLE_CHECK_PARAMS=$<BOOL:${BOYLE_CHECK_PARAMS}>)

set(CPM_SOURCE_CACHE ${CMAKE_SOURCE_DIR}/third_party CACHE PATH "CPM source cache")

if(BOYLE_BUILD_TESTING)
  include(CTest)
  enable_testing()
endif()

include(cmake/utils.cmake)
add_subdirectory(cmake/third_party)

add_subdirectory(src)
add_subdirectory(tests)

if(BOYLE_ENABLE_INSTALL)
  include(cmake/export.cmake)
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
  include(cmake/compile_commands.cmake)
endif()
