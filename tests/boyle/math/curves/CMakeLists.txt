boyle_cxx_test(
  NAME
    math_curve_proxy_test
  SRCS
    "curve_proxy_test.cpp"
  DEPS
    math_curve_proxy
    math_piecewise_linear_curve
    math_piecewise_cubic_curve
    math_piecewise_quintic_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_linear_curve2_test
  SRCS
    "piecewise_linear_curve2_test.cpp"
  DEPS
    math_piecewise_linear_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_cubic_curve2_test
  SRCS
    "piecewise_cubic_curve2_test.cpp"
  DEPS
    math_piecewise_cubic_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_quintic_curve2_test
  SRCS
    "piecewise_quintic_curve2_test.cpp"
  DEPS
    math_piecewise_quintic_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_linear_curve3_test
  SRCS
    "piecewise_linear_curve3_test.cpp"
  DEPS
    math_piecewise_linear_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_cubic_curve3_test
  SRCS
    "piecewise_cubic_curve3_test.cpp"
  DEPS
    math_piecewise_cubic_curve
)

boyle_cxx_test(
  NAME
    math_piecewise_quintic_curve3_test
  SRCS
    "piecewise_quintic_curve3_test.cpp"
  DEPS
    math_piecewise_quintic_curve
)
