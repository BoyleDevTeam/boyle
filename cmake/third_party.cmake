include(ExternalProject)

if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  find_package(LLVM 19.1.0 REQUIRED CONFIG)
endif()

find_package(Threads REQUIRED MODULE)

find_package(OpenMP REQUIRED MODULE)

include(${CMAKE_SOURCE_DIR}/cmake/OpenBLAS.cmake)

CPMAddPackage(
  NAME Boost
  VERSION 1.88.0
  URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.tar.xz
  URL_HASH SHA256=f48b48390380cfb94a629872346e3a81370dc498896f16019ade727ab72eb1ec
  OPTIONS
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "BUILD_SHARED_LIBS OFF"
    "BOOST_RUNTIME_LINK static"
    "BOOST_ENABLE_CMAKE ON"
    "BOOST_SKIP_INSTALL_RULES OFF"
    "BOOST_INCLUDE_LIBRARIES headers\\\;serialization"
)

CPMAddPackage(
  NAME cxxopts
  VERSION 3.2.1
  URL https://github.com/jarro2783/cxxopts/archive/refs/tags/v3.2.0.tar.gz
  URL_HASH SHA256=9f43fa972532e5df6c5fd5ad0f5bac606cdec541ccaf1732463d8070bbb7f03b
  OPTIONS
    "CXXOPTS_ENABLE_INSTALL ON"
)

CPMAddPackage(
  NAME fmt
  VERSION 11.2.0
  URL https://github.com/fmtlib/fmt/archive/refs/tags/11.2.0.tar.gz
  URL_HASH SHA256=bc23066d87ab3168f27cef3e97d545fa63314f5c79df5ea444d41d56f962c6af
  OPTIONS
    "BUILD_SHARED_LIBS OFF"
    "FMT_INSTALL ON"
)

CPMAddPackage(
  NAME spdlog
  VERSION 1.15.3
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.15.3.tar.gz
  URL_HASH SHA256=15a04e69c222eb6c01094b5c7ff8a249b36bb22788d72519646fb85feb267e67
  OPTIONS
    "CMAKE_POLICY_VERSION_MINIMUM 3.10"
    "SPDLOG_BUILD_SHARED OFF"
    "SPDLOG_INSTALL ON"
    "SPDLOG_FMT_EXTERNAL_HO ON"
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

CPMAddPackage(
  NAME msft_proxy
  VERSION 3.3.0
  URL https://github.com/microsoft/proxy/archive/refs/tags/3.3.0.tar.gz
  URL_HASH SHA256=9a5e89e70082cbdd937e80f5113f4ceb47bf6361cf7b88cb52782906a1b655cc
  OPTIONS
    "BUILD_TESTING OFF"
)

include(${CMAKE_SOURCE_DIR}/cmake/pocketfft.cmake)

CPMAddPackage(
  NAME qdldl
  VERSION 0.1.8
  URL https://github.com/osqp/qdldl/archive/refs/tags/v0.1.8.tar.gz
  URL_HASH SHA256=ecf113fd6ad8714f16289eb4d5f4d8b27842b6775b978c39def5913f983f6daa
  OPTIONS
    "CMAKE_POLICY_VERSION_MINIMUM 3.10"
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
    "OSQP_BUILD_SHARED_LIB OFF"
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
    "MATPLOTPP_BUILD_SHARED_LIBS OFF"
    "MATPLOTPP_BUILD_TESTS OFF"
    "MATPLOTPP_BUILD_EXAMPLES OFF"
    "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-ignored-attributes"
)
