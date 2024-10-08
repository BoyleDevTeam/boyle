if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  find_package(LLVM REQUIRED CONFIG)
endif()

CPMAddPackage(
  NAME cxxopts
  GITHUB_REPOSITORY "jarro2783/cxxopts"
  GIT_TAG "v3.2.1"
)

CPMAddPackage(
  NAME spdlog
  GITHUB_REPOSITORY "gabime/spdlog"
  GIT_TAG "v1.14.1"
  OPTIONS
    "SPDLOG_ENABLE_PCH ON"
    "SPDLOG_BUILD_PIC ON"
    "SPDLOG_TIDY OFF"
    "SPDLOG_USE_STD_FORMAT ON"
)

add_compile_definitions(SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)

CPMAddPackage(
  NAME doctest
  GITHUB_REPOSITORY "doctest/doctest"
  GIT_TAG "v2.4.9"
  OPTIONS
    "DOCTEST_WITH_MAIN_IN_STATIC_LIB ON"
    "DOCTEST_NO_INSTALL ON"
    "DOCTEST_USE_STD_HEADERS ON"
)

add_compile_definitions(DOCTEST_CONFIG_SUPER_FAST_ASSERTS)

# CPMAddPackage(
#   NAME nlohmann_json
#   GITHUB_REPOSITORY "nlohmann/json"
#   GIT_TAG "v3.11.2"
# )

CPMAddPackage(
  NAME msft_proxy
  GITHUB_REPOSITORY "microsoft/proxy"
  GIT_TAG "9636e15b359bba6ce588412f4db15661a39b8c47"
  OPTIONS
    "BUILD_TESTING OFF"
)

find_package(Boost 1.79.0 REQUIRED
  COMPONENTS
    serialization
)

# CPMAddPackage(
#   NAME OpenBLAS
#   GITHUB_REPOSITORY "OpenMathLib/OpenBLAS"
#   GIT_TAG "v0.3.27"
#   OPTIONS
#     "BUILD_STATIC_LIBS ON"
# )

# CPMAddPackage(
#   NAME autodiff
#   GITHUB_REPOSITORY "autodiff/autodiff"
#   GIT_TAG "v1.1.2"
#   OPTIONS
#     "AUTODIFF_BUILD_TESTS OFF"
#     "AUTODIFF_BUILD_PYTHON OFF"
#     "AUTODIFF_BUILD_EXAMPLES OFF"
#     "AUTODIFF_BUILD_DOCS OFF"
# )

CPMAddPackage(
  NAME osqp
  GITHUB_REPOSITORY "osqp/osqp"
  GIT_TAG "master"
  FORCE True
  OPTIONS
    "QDLDL_LONG OFF"
    "QDLDL_FLOAT OFF"
    "OSQP_USE_LONG OFF"
    "OSQP_USE_FLOAT OFF"
    "OSQP_BUILD_DEMO_EXE OFF"
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "CMAKE_POLICY_DEFAULT_CMP0169 OLD"
)

CPMAddPackage(
  NAME Matplot++
  GITHUB_REPOSITORY "alandefreitas/matplotplusplus"
  GIT_TAG "v1.2.1"
  OPTIONS
    "CMAKE_CXX_STANDARD 17"
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
)
