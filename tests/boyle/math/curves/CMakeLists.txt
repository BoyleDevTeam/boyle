add_subdirectory(piecewise_curves)

boyle_cxx_test(
  NAME
    math_curve2_proxy_test
  SRCS
    "curve2_proxy_test.cpp"
  DEPS
    math_curve2_proxy
    math_piecewise_linear_curve2
    math_piecewise_cubic_curve2
    math_piecewise_quintic_curve2
)
