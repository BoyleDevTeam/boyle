if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  find_package(LLVM REQUIRED CONFIG)
endif()

find_package(Threads REQUIRED MODULE)

find_package(OpenMP REQUIRED MODULE)

find_package(Boost 1.86.0 REQUIRED CONFIG
  COMPONENTS
    serialization
)

find_package(OpenBLAS REQUIRED CONFIG)

find_package(Eigen3 3.4.0 REQUIRED CONFIG)

set_target_properties(Eigen3::Eigen
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "EIGEN_USE_BLAS;EIGEN_USE_LAPACKE"
)

CPMAddPackage(
  NAME cxxopts
  GITHUB_REPOSITORY "jarro2783/cxxopts"
  GIT_TAG "v3.1.1"
  OPTIONS
    "CXXOPTS_ENABLE_INSTALL ON"
)

CPMAddPackage(
  NAME fmt
  GITHUB_REPOSITORY "fmtlib/fmt"
  GIT_TAG "9.1.0"
  OPTIONS
    "FMT_INSTALL ON"
)

CPMAddPackage(
  NAME spdlog
  GITHUB_REPOSITORY "gabime/spdlog"
  GIT_TAG "v1.12.0"
  OPTIONS
    "SPDLOG_INSTALL ON"
    "SPDLOG_FMT_EXTERNAL_HO ON"
)

if(spdlog_ADDED)
  set_target_properties(spdlog spdlog_header_only
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE"
  )
else()
  set_target_properties(spdlog::spdlog spdlog::spdlog_header_only
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE"
  )
endif()

CPMAddPackage(
  NAME doctest
  GITHUB_REPOSITORY "doctest/doctest"
  GIT_TAG "v2.4.11"
  OPTIONS
    "DOCTEST_WITH_MAIN_IN_STATIC_LIB ON"
    "DOCTEST_NO_INSTALL ON"
)

if(doctest_ADDED)
  set_target_properties(doctest
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "DOCTEST_CONFIG_SUPER_FAST_ASSERTS"
  )
else()
  set_target_properties(doctest::doctest
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "DOCTEST_CONFIG_SUPER_FAST_ASSERTS"
  )
endif()

CPMAddPackage(
  NAME msft_proxy
  GITHUB_REPOSITORY "microsoft/proxy"
  GIT_TAG "3.2.0"
  OPTIONS
    "BUILD_TESTING OFF"
)

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
  NAME pocketfft
  GITHUB_REPOSITORY "mreineck/pocketfft"
  GIT_TAG "cpp"
  DOWNLOAD_ONLY True
)

add_library(pocketfft INTERFACE)
target_sources(pocketfft
  INTERFACE
    FILE_SET HEADERS
      TYPE HEADERS
      BASE_DIRS ${pocketfft_SOURCE_DIR}
      FILES ${pocketfft_SOURCE_DIR}/pocketfft_hdronly.h
)
target_include_directories(pocketfft
  INTERFACE
    $<BUILD_INTERFACE:${pocketfft_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_definitions(pocketfft
  INTERFACE
    "POCKETFFT_NO_MULTITHREADING"
)

install(
  TARGETS
    pocketfft
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

CPMAddPackage(
  NAME osqp
  GITHUB_REPOSITORY "osqp/osqp"
  GIT_TAG "c7de4a6748e31be0b91a2cac2eb51625c89ca380"
  FORCE True
  OPTIONS
    "OSQP_USE_LONG OFF"
    "OSQP_USE_FLOAT OFF"
    "OSQP_BUILD_DEMO_EXE OFF"
    "OSQP_ENABLE_PROFILING OFF"
    "OSQP_ENABLE_INTERRUPT OFF"
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "CMAKE_POLICY_DEFAULT_CMP0069 NEW"
    "CMAKE_POLICY_DEFAULT_CMP0169 OLD"
)

CPMAddPackage(
  NAME Matplot++
  GITHUB_REPOSITORY "alandefreitas/matplotplusplus"
  GIT_TAG "v1.2.2"
  OPTIONS
    "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-ignored-attributes"
)
