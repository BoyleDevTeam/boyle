enable_language(Fortran)

set(_OPENBLAS_CC ${CMAKE_C_COMPILER_LAUNCHER}\ ${CMAKE_C_COMPILER})
set(_OPENBLAS_FC ${CMAKE_Fortran_COMPILER_LAUNCHER}\ ${CMAKE_Fortran_COMPILER})

execute_process(COMMAND xcrun --find ranlib OUTPUT_VARIABLE _OPENBLAS_RANLIB OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND xcrun --find ar OUTPUT_VARIABLE _OPENBLAS_AR OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" GCC_MAJOR_VERSION "${CMAKE_C_COMPILER_VERSION}")
  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "")
  file(GLOB_RECURSE _GFORTRAN_LIBRARY
    "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.a"
    "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.dylib"
  )
  if(_GFORTRAN_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
  endif()
  list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
  file(GLOB_RECURSE _QUADMATH_LIBRARY
    "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.a"
    "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.dylib"
  )
  if(_QUADMATH_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find quadmath library for OpenBLAS")
  endif()
  list(GET _QUADMATH_LIBRARY 0 _QUADMATH_LIBRARY)
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")
  unset(_GFORTRAN_LIBRARY)
  unset(_QUADMATH_LIBRARY)

elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")

  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "-rtlib=compiler-rt")
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "")
  file(GLOB _OPENBLAS_INTERFACE_LINK_LIBRARIES
    "/opt/homebrew/opt/flang/lib/clang/${LLVM_VERSION_MAJOR}/lib/darwin/libflang_rt.runtime.dylib"
    "/opt/homebrew/opt/flang/lib/libFortranDecimal.dylib"
  )
  if(_OPENBLAS_INTERFACE_LINK_LIBRARIES STREQUAL "")
    message(FATAL_ERROR "Could not find fortran libraries for OpenBLAS")
  endif()

elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" GCC_MAJOR_VERSION "${CMAKE_C_COMPILER_VERSION}")
  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "-rtlib=compiler-rt")
  file(GLOB_RECURSE _GFORTRAN_LIBRARY
    "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.a"
    "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.dylib"
  )
  if(_GFORTRAN_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find gfortran library for OpenBLAS")
  endif()
  list(GET _GFORTRAN_LIBRARY 0 _GFORTRAN_LIBRARY)
  file(GLOB_RECURSE _QUADMATH_LIBRARY
    "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.a"
    "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.dylib"
  )
  if(_QUADMATH_LIBRARY STREQUAL "")
    message(FATAL_ERROR "Could not find quadmath library for OpenBLAS")
  endif()
  list(GET _QUADMATH_LIBRARY 0 _QUADMATH_LIBRARY)
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")
  unset(_GFORTRAN_LIBRARY)
  unset(_QUADMATH_LIBRARY)

endif()

ExternalProject_Add(OpenBLAS
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/_deps
  TMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-tmp
  STAMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-stamp
  SOURCE_DIR ${OpenBLAS_SOURCE_DIR}
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build
  INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-install
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different <SOURCE_DIR>/. <BINARY_DIR> && ${CMAKE_COMMAND} -E env CC=${_OPENBLAS_CC} FC=${_OPENBLAS_FC} ASMFLAGS=${_OPENBLAS_ASMFLAGS} CFLAGS=${_OPENBLAS_CFLAGS} FFLAGS=${_OPENBLAS_FFLAGS} LDFLAGS=${_OPENBLAS_LDFLAGS} DYNAMIC_ARCH=0 NO_AFFINITY=0 USE_THREAD=0 USE_LOCKING=0 USE_OPENMP=0 make RANLIB=${_OPENBLAS_RANLIB} AR=${_OPENBLAS_AR} > /dev/null 2>&1
  INSTALL_COMMAND make -C <BINARY_DIR> PREFIX=<INSTALL_DIR> install > /dev/null 2>&1
  BUILD_BYPRODUCTS
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/libopenblas.dylib
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-build/libopenblas.a
  INSTALL_BYPRODUCTS
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-install/lib/libopenblas.dylib
    ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-install/lib/libopenblas.a
)

set(_OPENBLAS_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/openblas-install)
set(_OPENBLAS_INCLUDE_DIR ${_OPENBLAS_INSTALL_DIR}/include)
set(_OPENBLAS_SHARED_LIB ${_OPENBLAS_INSTALL_DIR}/lib/libopenblas.dylib)
set(_OPENBLAS_STATIC_LIB ${_OPENBLAS_INSTALL_DIR}/lib/libopenblas.a)

unset(_OPENBLAS_CC)
unset(_OPENBLAS_FC)
unset(_OPENBLAS_RANLIB)
unset(_OPENBLAS_AR)
unset(_OPENBLAS_ASMFLAGS)
unset(_OPENBLAS_CFLAGS)
unset(_OPENBLAS_FFLAGS)
unset(_OPENBLAS_LDFLAGS)

add_library(blas_lapack INTERFACE)
add_library(blas_lapack_static INTERFACE)

add_dependencies(blas_lapack OpenBLAS)
add_dependencies(blas_lapack_static OpenBLAS)

set_target_properties(blas_lapack
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK"
    INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${_OPENBLAS_INCLUDE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    INTERFACE_LINK_LIBRARIES "$<BUILD_INTERFACE:${_OPENBLAS_SHARED_LIB}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.dylib>;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
)
set_target_properties(blas_lapack_static
  PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BLAS_LAPACK"
    INTERFACE_COMPILE_OPTIONS "-Wno-c99-extensions"
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${_OPENBLAS_INCLUDE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    INTERFACE_LINK_LIBRARIES "$<BUILD_INTERFACE:${_OPENBLAS_STATIC_LIB}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.a>;${_OPENBLAS_INTERFACE_LINK_LIBRARIES}"
)

if(BOYLE_ENABLE_INSTALL)
  install(
    DIRECTORY ${_OPENBLAS_INSTALL_DIR}/.
    DESTINATION ${CMAKE_INSTALL_PREFIX}
  )
  install(
    TARGETS blas_lapack blas_lapack_static
    EXPORT ${CMAKE_PROJECT_NAME}Targets
  )
endif()

set(BLAS_FOUND True CACHE INTERNAL "set BLAS found")
set(BLAS_INCLUDE_DIRS "$<BUILD_INTERFACE:${_OPENBLAS_INCLUDE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>" CACHE INTERNAL "set BLAS include directory")
set(BLAS_LIBRARIES "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>" CACHE INTERNAL "set BLAS library")
set(BLAS_LINKER_FLAGS "" CACHE INTERNAL "set BLAS linker flags")
set(LAPACK_FOUND True CACHE INTERNAL "set LAPACK found")
set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS} CACHE INTERNAL "set LAPACK include directory")
set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "set LAPACK library")
set(LAPACK_LINKER_FLAGS "${BLAS_LINKER_FLAGS}" CACHE INTERNAL "set LAPACK linker flags")

unset(_OPENBLAS_INSTALL_DIR)
unset(_OPENBLAS_INCLUDE_DIR)
unset(_OPENBLAS_SHARED_LIB)
unset(_OPENBLAS_STATIC_LIB)
unset(_OPENBLAS_INTERFACE_LINK_LIBRARIES)
