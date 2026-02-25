add_library(blas_lapack INTERFACE)
add_library(blas_lapack_static INTERFACE)
if(BOYLE_USE_BLAS_LAPACK STREQUAL "OpenBLAS")
  CPMAddPackage(
    NAME OpenBLAS
    VERSION 0.3.33
    URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.33/OpenBLAS-0.3.33.tar.gz
    URL_HASH SHA256=6761af1d9f5d353ab4f0b7497be2643313b36c8f31caec0144bfef198e71e6ab
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
