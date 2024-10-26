add_subdirectory(curves)
add_subdirectory(functions)
add_subdirectory(mdfunctions)
add_subdirectory(sparse_matrix)

boyle_cxx_test(
  NAME
    math_eigen_serialization_test
  SRCS
    "eigen_serialization_test.cpp"
  DEPS
    math_eigen_serialization
)

boyle_cxx_test(
  NAME
    math_fft_test
  SRCS
    "fft_test.cpp"
  DEPS
    math_fft
)

boyle_cxx_test(
  NAME
    math_chebyshev_test
  SRCS
    "chebyshev_test.cpp"
  DEPS
    math_chebyshev
    math_utils
)

boyle_cxx_test(
  NAME
    math_vec2_test
  SRCS
    "vec2_test.cpp"
  DEPS
    fmt::fmt-header-only
    math_vec2
)

boyle_cxx_test(
  NAME
    math_vec3_test
  SRCS
    "vec3_test.cpp"
  DEPS
    fmt::fmt-header-only
    math_vec3
)

boyle_cxx_test(
  NAME
    math_cubic_interpolation_test
  SRCS
    "cubic_interpolation_test.cpp"
  DEPS
    math_cubic_interpolation
)

boyle_cxx_test(
  NAME
    math_quintic_interpolation_test
  SRCS
    "quintic_interpolation_test.cpp"
  DEPS
    math_quintic_interpolation
)
