enable_language(Fortran)

string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" _FORTRAN_MAJOR_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")

set(_OPENBLAS_CC ${CMAKE_C_COMPILER_LAUNCHER}\ ${CMAKE_C_COMPILER})
set(_OPENBLAS_FC ${CMAKE_Fortran_COMPILER_LAUNCHER}\ ${CMAKE_Fortran_COMPILER})

execute_process(COMMAND xcrun --find ranlib OUTPUT_VARIABLE _OPENBLAS_RANLIB OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND xcrun --find ar OUTPUT_VARIABLE _OPENBLAS_AR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(_GFORTRAN_GLOBS
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.a"
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.dylib"
)
set(_QUADMATH_GLOBS
  "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.a"
  "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.dylib"
)

if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${_OPENBLAS_ASMFLAGS})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "")
  boyle_glob_first_lib(_GFORTRAN_LIBRARY ${_GFORTRAN_GLOBS})
  boyle_glob_first_lib(_QUADMATH_LIBRARY ${_QUADMATH_GLOBS})
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")

elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")

  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${_OPENBLAS_ASMFLAGS})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "-rtlib=compiler-rt")
  file(GLOB _OPENBLAS_INTERFACE_LINK_LIBRARIES
    "/opt/homebrew/opt/flang/lib/clang/${_FORTRAN_MAJOR_VERSION}/lib/darwin/libflang_rt.runtime.dylib"
    "/opt/homebrew/opt/flang/lib/libFortranDecimal.dylib"
  )
  if(_OPENBLAS_INTERFACE_LINK_LIBRARIES STREQUAL "")
    message(FATAL_ERROR "Could not find fortran libraries for OpenBLAS")
  endif()

elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  set(_OPENBLAS_ASMFLAGS ${CMAKE_C_FLAGS}\ ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_CFLAGS ${_OPENBLAS_ASMFLAGS})
  set(_OPENBLAS_FFLAGS ${CMAKE_Fortran_FLAGS}\ ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}})
  set(_OPENBLAS_LDFLAGS "-rtlib=compiler-rt")
  boyle_glob_first_lib(_GFORTRAN_LIBRARY ${_GFORTRAN_GLOBS})
  boyle_glob_first_lib(_QUADMATH_LIBRARY ${_QUADMATH_GLOBS})
  set(_OPENBLAS_INTERFACE_LINK_LIBRARIES "${_GFORTRAN_LIBRARY};${_QUADMATH_LIBRARY}")

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

add_dependencies(blas_lapack OpenBLAS)
add_dependencies(blas_lapack_static OpenBLAS)

target_compile_definitions(blas_lapack INTERFACE BOYLE_USE_BLAS_LAPACK)
target_compile_options(blas_lapack INTERFACE -Wno-c99-extensions)
target_include_directories(blas_lapack INTERFACE
  "$<BUILD_INTERFACE:${_OPENBLAS_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(blas_lapack INTERFACE
  "$<BUILD_INTERFACE:${_OPENBLAS_SHARED_LIB}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.dylib>"
  ${_OPENBLAS_INTERFACE_LINK_LIBRARIES}
)

target_compile_definitions(blas_lapack_static INTERFACE BOYLE_USE_BLAS_LAPACK)
target_compile_options(blas_lapack_static INTERFACE -Wno-c99-extensions)
target_include_directories(blas_lapack_static INTERFACE
  "$<BUILD_INTERFACE:${_OPENBLAS_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(blas_lapack_static INTERFACE
  "$<BUILD_INTERFACE:${_OPENBLAS_STATIC_LIB}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/libopenblas.a>"
  ${_OPENBLAS_INTERFACE_LINK_LIBRARIES}
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
