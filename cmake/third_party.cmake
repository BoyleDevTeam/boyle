set(CMAKE_COMPILE_WARNING_AS_ERROR False)

CPMAddPackage(
  NAME cxxopts
  GITHUB_REPOSITORY "jarro2783/cxxopts"
  GIT_TAG "v3.1.1"
  OPTIONS
    "CXXOPTS_BUILD_EXAMPLES NO"
    "CXXOPTS_BUILD_TESTS NO"
    "CXXOPTS_ENABLE_INSTALL YES"
    "CXXOPTS_USE_UNICODE_HELP YES"
)

CPMAddPackage(
  NAME spdlog
  GITHUB_REPOSITORY "gabime/spdlog"
  GIT_TAG "v1.12.0"
  OPTIONS
    "SPDLOG_ENABLE_PCH ON"
    "SPDLOG_BUILD_PIC ON"
    "SPDLOG_TIDY OFF"
    "SPDLOG_USE_STD_FORMAT ON"
)

CPMAddPackage(
  NAME doctest
  GITHUB_REPOSITORY "doctest/doctest"
  GIT_TAG "v2.4.9"
)

CPMAddPackage(
  NAME nlohmann_json
  GITHUB_REPOSITORY "nlohmann/json"
  GIT_TAG "v3.11.2"
  OPTIONS
    "JSON_BuildTests OFF"
)

# CPMAddPackage(
#   NAME Taskflow
#   GITHUB_REPOSITORY "taskflow/taskflow"
#   GIT_TAG "v3.6.0"
#   OPTIONS
#     "TF_BUILD_TESTS OFF"
#     "TF_BUILD_EXAMPLES OFF"
# )

# CPMAddPackage(
#  NAME zpp_bits
#  GITHUB_REPOSITORY "eyalz800/zpp_bits"
#  GIT_TAG "v4.4.17"
#  FORCE True
#  DOWNLOAD_ONLY True
# )

# if(zpp_bits_ADDED)
#  add_library(zpp_bits INTERFACE IMPORTED)
#  target_sources(zpp_bits INTERFACE "${zpp_bits_SOURCE_DIR}/zpp_bits.h")
#  target_include_directories(zpp_bits INTERFACE "${zpp_bits_SOURCE_DIR}")
#  target_compile_features(zpp_bits INTERFACE cxx_std_20)
# endif()

find_package(Boost
  REQUIRED
  COMPONENTS
    serialization
)

CPMAddPackage(
  NAME Eigen3
  GITLAB_REPOSITORY "libeigen/eigen"
  GIT_TAG "3.4.0"
  OPTIONS
    "BUILD_TESTING OFF"
    "EIGEN_TEST_CXX11 OFF"
)

CPMAddPackage(
  NAME osqp
  GITHUB_REPOSITORY "osqp/osqp"
  GIT_TAG "v1.0.0.beta1"
  FORCE True
  OPTIONS
    "QDLDL_LONG OFF"
    "QDLDL_FLOAT OFF"
    "OSQP_USE_LONG OFF"
    "OSQP_USE_FLOAT OFF"
    "OSQP_BUILD_SHARED_LIB OFF"
    "OSQP_BUILD_DEMO_EXE OFF"
)

# CPMAddPackage(
#   NAME foonathan_memory
#   GITHUB_REPOSITORY "foonathan/memory"
#   GIT_TAG "v0.7-3"
#   OPTIONS
#     "FOONATHAN_MEMORY_BUILD_EXAMPLES OFF"
#     "FOONATHAN_MEMORY_BUILD_TESTS OFF"
#     "BUILD_SHARED_LIBS OFF"
# )

# CPMAddPackage(
#   NAME fastcdr
#   GITHUB_REPOSITORY "eProsima/Fast-CDR"
#   GIT_TAG "v1.1.0"
#   OPTIONS
#     "EPROSIMA_BUILD_TESTS OFF"
#     "BUILD_SHARED_LIBS OFF"
#     "FORCE_CXX 17"
# )

# CPMAddPackage(
#   NAME fastrtps
#   GITHUB_REPOSITORY "eProsima/Fast-DDS"
#   GIT_TAG "v2.11.2"
#   OPTIONS
#     "EPROSIMA_BUILD_TESTS OFF"
#     "BUILD_SHARED_LIBS OFF"
#     "FORCE_CXX 17"
#     "THIRDPARTY ON"
# )

CPMAddPackage(
  NAME Matplot++
  GITHUB_REPOSITORY "alandefreitas/matplotplusplus"
  GIT_TAG "v1.2.0"
  OPTIONS
    "MATPLOTPP_BUILD_EXAMPLES OFF"
    "MATPLOTPP_BUILD_TESTS OFF"
    "MATPLOTPP_BUILD_SHARED_LIBS OFF"
    "CMAKE_INTERPROCEDURAL_OPTIMIZATION ON"
    "CMAKE_CXX_STANDARD 17"
)

set(CMAKE_COMPILE_WARNING_AS_ERROR True)
