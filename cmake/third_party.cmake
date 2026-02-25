include(ExternalProject)

if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  find_package(LLVM 21.1.6 REQUIRED CONFIG)
endif()

find_package(Threads REQUIRED MODULE)

find_package(OpenMP REQUIRED MODULE)

CPMAddPackage(
  NAME Boost
  VERSION 1.90.0
  URL https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-cmake.tar.xz
  URL_HASH SHA256=aca59f889f0f32028ad88ba6764582b63c916ce5f77b31289ad19421a96c555f
  OPTIONS
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "CMAKE_UNITY_BUILD OFF"
    "BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS}"
    "BOOST_RUNTIME_LINK $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
    "BOOST_ENABLE_CMAKE ON"
    "BOOST_SKIP_INSTALL_RULES OFF"
    "BOOST_INCLUDE_LIBRARIES headers\\\;serialization"
)

CPMAddPackage(
  NAME cxxopts
  VERSION 3.3.1
  URL https://github.com/jarro2783/cxxopts/archive/refs/tags/v3.3.1.tar.gz
  URL_HASH SHA256=3bfc70542c521d4b55a46429d808178916a579b28d048bd8c727ee76c39e2072
  OPTIONS
    "CXXOPTS_ENABLE_INSTALL ON"
)

CPMAddPackage(
  NAME spdlog
  VERSION 1.17.0
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.17.0.tar.gz
  URL_HASH SHA256=d8862955c6d74e5846b3f580b1605d2428b11d97a410d86e2fb13e857cd3a744
  OPTIONS
    "SPDLOG_BUILD_SHARED ${BUILD_SHARED_LIBS}"
    "SPDLOG_INSTALL ON"
    "SPDLOG_USE_STD_FORMAT ON"
)
set_target_properties(spdlog spdlog_header_only
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE"
)

CPMAddPackage(
  NAME doctest
  VERSION 2.4.12
  URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.12.tar.gz
  URL_HASH SHA256=73381c7aa4dee704bd935609668cf41880ea7f19fa0504a200e13b74999c2d70
  PATCHES
    ${CMAKE_SOURCE_DIR}/cmake/patches/doctest.patch
  OPTIONS
    "CMAKE_POLICY_VERSION_MINIMUM 3.10"
    "CMAKE_POLICY_DEFAULT_CMP0069 NEW"
    "DOCTEST_WITH_MAIN_IN_STATIC_LIB ON"
    "DOCTEST_NO_INSTALL ON"
)
set_target_properties(doctest
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "DOCTEST_CONFIG_SUPER_FAST_ASSERTS"
)

include(${CMAKE_SOURCE_DIR}/cmake/msft_proxy.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/pocketfft.cmake)

if(BOYLE_USE_BLAS_LAPACK STREQUAL "OpenBLAS")
  include(${CMAKE_SOURCE_DIR}/cmake/OpenBLAS.cmake)
elseif(BOYLE_USE_BLAS_LAPACK STREQUAL "OFF")
  add_library(blas_lapack INTERFACE)
  add_library(blas_lapack_static INTERFACE)
endif()

CPMAddPackage(
  NAME qdldl
  VERSION 0.1.9
  URL https://github.com/osqp/qdldl/archive/refs/tags/v0.1.9.tar.gz
  URL_HASH SHA256=7d1285b2db15cf2730dc83b3d16ed28412f558591108cca4f28d4438bf72ceed
  OPTIONS
    "CMAKE_POLICY_DEFAULT_CMP0069 NEW"
    "QDLDL_BUILD_STATIC_LIB OFF"
    "QDLDL_BUILD_SHARED_LIB OFF"
    "QDLDL_LONG OFF"
    "QDLDL_FLOAT OFF"
    "QDLDL_BUILD_DEMO_EXE OFF"
    "QDLDL_UNITTESTS OFF"
)

CPMAddPackage(
  NAME osqp
  VERSION 1.0.0
  URL https://github.com/osqp/osqp/releases/download/v1.0.0/osqp-v1.0.0-src.tar.gz
  URL_HASH SHA256=ec0bb8fd34625d0ea44274ab3e991aa56e3e360ba30935ae62476557b101c646
  OPTIONS
    "CMAKE_POLICY_VERSION_MINIMUM 3.10"
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "OSQP_BUILD_SHARED_LIB ${BUILD_SHARED_LIBS}"
    "OSQP_USE_LONG OFF"
    "OSQP_USE_FLOAT OFF"
    "OSQP_BUILD_DEMO_EXE OFF"
    "OSQP_ENABLE_PRINTING OFF"
    "OSQP_ENABLE_PROFILING OFF"
    "OSQP_ENABLE_INTERRUPT OFF"
)

CPMAddPackage(
  NAME Matplot++
  VERSION 1.2.2
  URL https://github.com/alandefreitas/matplotplusplus/archive/refs/tags/v1.2.2.tar.gz
  URL_HASH SHA256=c7434b4fea0d0cc3508fd7104fafbb2fa7c824b1d2ccc51c52eaee26fc55a9a0
  PATCHES
    ${CMAKE_SOURCE_DIR}/cmake/patches/matplot++.patch
  OPTIONS
    "CMAKE_CXX_STANDARD 23"
    "MATPLOTPP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS}"
    "MATPLOTPP_BUILD_TESTS OFF"
    "MATPLOTPP_BUILD_EXAMPLES OFF"
)
