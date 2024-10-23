if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  find_package(LLVM REQUIRED CONFIG)
endif()

find_package(Threads REQUIRED)

find_package(Boost 1.79.0 REQUIRED
  COMPONENTS
    serialization
)

CPMAddPackage(
  NAME cxxopts
  GITHUB_REPOSITORY "jarro2783/cxxopts"
  GIT_TAG "v3.2.1"
  OPTIONS
    "CXXOPTS_ENABLE_INSTALL ON"
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
    "SPDLOG_INSTALL ON"
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
  GIT_TAG "v2.4.9"
  OPTIONS
    "DOCTEST_WITH_MAIN_IN_STATIC_LIB ON"
    "DOCTEST_NO_INSTALL ON"
    "DOCTEST_USE_STD_HEADERS ON"
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
  NAME nlohmann_json
  GITHUB_REPOSITORY "nlohmann/json"
  GIT_TAG "v3.11.2"
  OPTIONS
    "JSON_Install ON"
)

CPMAddPackage(
  NAME msft_proxy
  GITHUB_REPOSITORY "microsoft/proxy"
  GIT_TAG "3.0.0"
  OPTIONS
    "BUILD_TESTING OFF"
)

# CPMAddPackage(
#   NAME OpenBLAS
#   GITHUB_REPOSITORY "OpenMathLib/OpenBLAS"
#   GIT_TAG "v0.3.27"
#   OPTIONS
#     "BUILD_STATIC_LIBS ON"
# )

CPMAddPackage(
  NAME Eigen3
  GITLAB_REPOSITORY "libeigen/eigen"
  GIT_TAG "3.4"
  OPTIONS
    "EIGEN_BUILD_TESTING OFF"
    "EIGEN_BUILD_DOC OFF"
    "EIGEN_BUILD_DEMOS OFF"
    "EIGEN_BUILD_BLAS ON"
    "EIGEN_BUILD_LAPACK ON"
    "EIGEN_BUILD_CMAKE_PACKAGE ON"
)

if(Eigen3_ADDED)
  set_target_properties(eigen
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "EIGEN_USE_BLAS;EIGEN_USE_LAPACK"
  )
else()
  set_target_properties(Eigen3::Eigen
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "EIGEN_USE_BLAS;EIGEN_USE_LAPACK"
  )
endif()

CPMAddPackage(
  NAME autodiff
  GITHUB_REPOSITORY "autodiff/autodiff"
  GIT_TAG "v1.1.2"
  OPTIONS
    "AUTODIFF_BUILD_TESTS OFF"
    "AUTODIFF_BUILD_PYTHON OFF"
    "AUTODIFF_BUILD_EXAMPLES OFF"
    "AUTODIFF_BUILD_DOCS OFF"
)

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

install(
  TARGETS
    pocketfft
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

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
    "CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -Wno-cpp -Wno-pedantic"
    "CMAKE_POLICY_DEFAULT_CMP0069 NEW"
    "CMAKE_POLICY_DEFAULT_CMP0169 OLD"
)

CPMAddPackage(
  NAME Matplot++
  GITHUB_REPOSITORY "alandefreitas/matplotplusplus"
  GIT_TAG "v1.2.1"
  OPTIONS
    "CMAKE_CXX_STANDARD 17"
    "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-uninitialized -Wno-ignored-attributes"
)
