boyle_cxx_test(
  NAME
    math_piecewise_linear_function1_test
  SRCS
    "piecewise_linear_function1_test.cpp"
  DEPS
    math_piecewise_linear_function1
)

boyle_cxx_test(
  NAME
    math_piecewise_cubic_function1_test
  SRCS
    "piecewise_cubic_function1_test.cpp"
  DEPS
    math_piecewise_cubic_function1
)

boyle_cxx_test(
  NAME
    math_piecewise_quintic_function1_test
  SRCS
    "piecewise_quintic_function1_test.cpp"
  DEPS
    math_piecewise_quintic_function1
)
