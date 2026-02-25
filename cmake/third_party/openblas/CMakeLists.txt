add_library(blas_lapack INTERFACE)
add_library(blas_lapack_static INTERFACE)
if(BOYLE_USE_BLAS_LAPACK STREQUAL "OpenBLAS")
  CPMAddPackage(
    NAME OpenBLAS
    VERSION 0.3.32
    GITHUB_REPOSITORY OpenMathLib/OpenBLAS
    GIT_TAG v0.3.32
    DOWNLOAD_ONLY True
  )
  if (WIN32)
    include(${CMAKE_CURRENT_SOURCE_DIR}/windows.cmake)
  elseif(APPLE)
    include(${CMAKE_CURRENT_SOURCE_DIR}/darwin.cmake)
  else()
    include(${CMAKE_CURRENT_SOURCE_DIR}/linux.cmake)
  endif()
endif()
