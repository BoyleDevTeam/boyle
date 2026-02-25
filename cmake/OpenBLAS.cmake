enable_language(Fortran)

set(_OPENBLAS_CC ${CMAKE_C_COMPILER_LAUNCHER}\ ${CMAKE_C_COMPILER})
set(_OPENBLAS_FC ${CMAKE_Fortran_COMPILER_LAUNCHER}\ ${CMAKE_Fortran_COMPILER})

if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" GCC_MAJOR_VERSION "${CMAKE_C_COMPILER_VERSION}")
  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
  set(_OPENBLAS_LDFLAGS -fuse-ld=bfd)
  file(GLOB_RECURSE _GFORTRAN_LIBRARY
    "/usr/lib/libgfortran.a"
    "/usr/lib/libgfortran.so"
    "/usr/local/lib64/libgfortran.a"
    "/usr/local/lib64/libgfortran.so"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libgfortran.a"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libgfortran.so"
  )
  if(_GFORTRAN_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
  endif()
  list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
  file(GLOB_RECURSE _QUADMATH_LIBRARY
    "/usr/lib/libquadmath.a"
    "/usr/lib/libquadmath.so"
    "/usr/local/lib64/libquadmath.a"
    "/usr/local/lib64/libquadmath.so"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libquadmath.a"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libquadmath.so"
  )
  if(_QUADMATH_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find quadmath library for OpenBLAS")
  endif()
  list(GET _QUADMATH_LIBRARY 0 _QUADMATH_LIBRARY)
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")
  unset(_GFORTRAN_LIBRARY)
  unset(_QUADMATH_LIBRARY)
elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=thin)
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=thin)
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=thin)
  set(_OPENBLAS_LDFLAGS -fuse-ld=lld\ -rtlib=compiler-rt)
  file(GLOB _OPENBLAS_INTERFACE_LINK_LIBRARIES
    "${LLVM_LIBRARY_DIR}/clang/${LLVM_VERSION_MAJOR}/lib/x86_64-pc-linux-gnu/libflang_rt.runtime.a"
    "${LLVM_LIBRARY_DIR}/libFortranDecimal.a"
  )
  if(_OPENBLAS_INTERFACE_LINK_LIBRARIES STREQUAL "")
    message(FATAL_ERROR "Could not find fortran libraries for OpenBLAS")
  endif()
elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" GCC_MAJOR_VERSION "${CMAKE_C_COMPILER_VERSION}")
  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=thin)
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=thin)
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\ -fno-lto)
  set(_OPENBLAS_LDFLAGS -fuse-ld=lld)
  file(GLOB_RECURSE _GFORTRAN_LIBRARY
    "/usr/lib/libgfortran.a"
    "/usr/lib/libgfortran.so"
    "/usr/local/lib64/libgfortran.a"
    "/usr/local/lib64/libgfortran.so"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libgfortran.a"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libgfortran.so"
  )
  if(_GFORTRAN_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
  endif()
  list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
  file(GLOB_RECURSE _QUADMATH_LIBRARY
    "/usr/lib/libquadmath.a"
    "/usr/lib/libquadmath.so"
    "/usr/local/lib64/libquadmath.a"
    "/usr/local/lib64/libquadmath.so"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libquadmath.a"
    "/usr/lib/gcc/x86_64-linux-gnu/${GCC_MAJOR_VERSION}/libquadmath.so"
  )
  if(_QUADMATH_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find quadmath library for OpenBLAS")
  endif()
  list(GET _QUADMATH_LIBRARY 0 _QUADMATH_LIBRARY)
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")
  unset(_GFORTRAN_LIBRARY)
  unset(_QUADMATH_LIBRARY)
endif()

CPMAddPackage(
  NAME OpenBLAS
  VERSION 0.3.31
  URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31.tar.gz
  URL_HASH SHA256=6dd2a63ac9d32643b7cc636eab57bf4e57d0ed1fff926dfbc5d3d97f2d2be3a6
  PATCHES
    ${CMAKE_SOURCE_DIR}/cmake/patches/openblas.patch
  DOWNLOAD_ONLY True
)

ExternalProject_Add(OpenBLAS
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/_deps
  TMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-tmp
  STAMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-stamp
  SOURCE_DIR ${OpenBLAS_SOURCE_DIR}
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different <SOURCE_DIR>/. <BINARY_DIR> && ${CMAKE_COMMAND} -E env CC=${_OPENBLAS_CC} FC=${_OPENBLAS_FC} ASMFLAGS=${_OPENBLAS_ASMFLAGS} CFLAGS=${_OPENBLAS_CFLAGS} FFLAGS=${_OPENBLAS_FFLAGS} LDFLAGS=${_OPENBLAS_LDFLAGS} DYNAMIC_ARCH=0 NO_AFFINITY=0 USE_THREAD=0 USE_LOCKING=0 USE_OPENMP=0 make > /dev/null 2>&1
  INSTALL_COMMAND make -C <BINARY_DIR> PREFIX=<INSTALL_DIR> install > /dev/null 2>&1
  BUILD_BYPRODUCTS
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/libopenblas.so
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/libopenblas.a
  INSTALL_BYPRODUCTS
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib/libopenblas.so
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib/libopenblas.a
)

add_library(blas_lapack INTERFACE)
add_library(blas_lapack_static INTERFACE)

add_dependencies(blas_lapack OpenBLAS)
add_dependencies(blas_lapack_static OpenBLAS)

set_target_properties(blas_lapack
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK"
    INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    INTERFACE_LINK_LIBRARIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib/libopenblas.so>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.so>;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
)
set_target_properties(blas_lapack_static
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK"
    INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    INTERFACE_LINK_LIBRARIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib/libopenblas.a>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.a>;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
)

if(BOYLE_ENABLE_INSTALL)
  install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/.
    DESTINATION ${CMAKE_INSTALL_PREFIX}
  )
  install(
    TARGETS blas_lapack blas_lapack_static
    EXPORT ${CMAKE_PROJECT_NAME}Targets
  )
endif()

set(BLAS_FOUND True CACHE INTERNAL "set BLAS found")
set(BLAS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>" CACHE INTERNAL "set BLAS include directory")
set(BLAS_LIBRARIES "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>" CACHE INTERNAL "set BLAS library")
set(BLAS_LINKER_FLAGS "" CACHE INTERNAL "set BLAS linker flags")
set(LAPACK_FOUND True CACHE INTERNAL "set LAPACK found")
set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS} CACHE INTERNAL "set LAPACK include directory")
set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "set LAPACK library")
set(LAPACK_LINKER_FLAGS "${BLAS_LINKER_FLAGS}" CACHE INTERNAL "set LAPACK linker flags")

unset(_OPENBLAS_CC)
unset(_OPENBLAS_FC)
unset(_OPENBLAS_ASMFLAGS)
unset(_OPENBLAS_CFLAGS)
unset(_OPENBLAS_FFLAGS)
unset(_OPENBLAS_LDFLAGS)
unset(_OPENBLAS_INTERFACE_LINK_LIBRARIES)
