add_library(openblas INTERFACE)
add_library(openblas_static INTERFACE)

if(BOYLE_USE_BLAS_LAPACK)
  set(_OPENBLAS_CC ${CMAKE_C_COMPILER_LAUNCHER}\ ${CMAKE_C_COMPILER})
  set(_OPENBLAS_FC ${CMAKE_Fortran_COMPILER_LAUNCHER}\ ${CMAKE_Fortran_COMPILER})

  if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
    set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
    set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\ -flto=auto\ -fno-fat-lto-objects)
    set(_OPENBLAS_LDFLAGS -fuse-ld=gold)
    file(GLOB_RECURSE _GFORTRAN_LIBRARY
      "/usr/lib*/libgfortran.a"
      "/usr/lib*/libgfortran.so.[0-9]"
    )
    if(_GFORTRAN_LIBRARY STREQUAL "")
      message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
    endif()
    list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
    file(GLOB_RECURSE _QUADMATH_LIBRARY
      "/usr/lib*/libquadmath.a"
      "/usr/lib*/libquadmath.so.[0-9]"
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
      "${LLVM_LIBRARY_DIR}/libFortran*.a"
    )
    if(_OPENBLAS_INTERFACE_LINK_LIBRARIES STREQUAL "")
      message(FATAL_ERROR "Could not find fortran libraries for OpenBLAS")
    endif()
  elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -fno-lto)
    set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}\ -fno-lto)
    set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}\ -fno-lto)
    set(_OPENBLAS_LDFLAGS -fuse-ld=gold)
    file(GLOB_RECURSE _GFORTRAN_LIBRARY
      "/usr/lib*/libgfortran.a"
      "/usr/lib*/libgfortran.so.[0-9]"
    )
    if(_GFORTRAN_LIBRARY STREQUAL "")
      message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
    endif()
    list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
    file(GLOB_RECURSE _QUADMATH_LIBRARY
      "/usr/lib*/libquadmath.a"
      "/usr/lib*/libquadmath.so.[0-9]"
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
    VERSION 0.3.30
    URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz
    URL_HASH SHA256=27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d
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
    BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different <SOURCE_DIR>/. <BINARY_DIR> && ${CMAKE_COMMAND} -E env CC=${_OPENBLAS_CC} FC=${_OPENBLAS_FC} ASMFLAGS=${_OPENBLAS_ASMFLAGS} CFLAGS=${_OPENBLAS_CFLAGS} FFLAGS=${_OPENBLAS_FFLAGS} LDFLAGS=${_OPENBLAS_LDFLAGS} ARCH=x86-64 DYNAMIC_ARCH=0 NO_AFFINITY=0 USE_THREAD=0 USE_LOCKING=0 USE_OPENMP=0 BUILD_LAPACK_DEPRECATED=1 make > /dev/null 2>&1
    INSTALL_COMMAND make -C <BINARY_DIR> PREFIX=<INSTALL_DIR> install > /dev/null 2>&1 && ${CMAKE_COMMAND} -E copy_directory_if_different <INSTALL_DIR>/. ${CMAKE_INSTALL_PREFIX}
  )

  add_dependencies(openblas OpenBLAS)
  add_dependencies(openblas_static OpenBLAS)

  set_target_properties(openblas
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK=1"
      INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
      INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
      INTERFACE_LINK_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>"
      INTERFACE_LINK_LIBRARIES "libopenblas.so;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
  )
  set_target_properties(openblas_static
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK=1"
      INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
      INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
      INTERFACE_LINK_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/lib>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>"
      INTERFACE_LINK_LIBRARIES "libopenblas.a;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
  )

  set(BLAS_FOUND True CACHE INTERNAL "set BLAS found")
  set(BLAS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/usr/include CACHE INTERNAL "set BLAS include directory")
  set(BLAS_LIBRARIES openblas_static CACHE INTERNAL "set BLAS library")
  set(BLAS_LINKER_FLAGS "${_OPENBLAS_INTERFACE_LINK_LIBRARIES}" CACHE INTERNAL "set BLAS linker flags")
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
endif()

if(BOYLE_ENABLE_INSTALL)
  install(
    TARGETS
      openblas openblas_static
    EXPORT ${CMAKE_PROJECT_NAME}Targets
  )
endif()
