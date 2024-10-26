boyle_cxx_test(
  NAME
    math_piecewise_linear_curve2_test
  SRCS
    "piecewise_linear_curve2_test.cpp"
  DEPS
    math_piecewise_linear_curve2
)

boyle_cxx_test(
  NAME
    math_piecewise_cubic_curve2_test
  SRCS
    "piecewise_cubic_curve2_test.cpp"
  DEPS
    math_piecewise_cubic_curve2
)

boyle_cxx_test(
  NAME
    math_piecewise_quintic_curve2_test
  SRCS
    "piecewise_quintic_curve2_test.cpp"
  DEPS
    math_piecewise_quintic_curve2
)
