add_boyle_test(math_curve2_proxy curve2_proxy_test.cpp)
target_link_libraries(math_curve2_proxy_test
  PRIVATE
    math_piecewise_linear_curve2
    math_piecewise_cubic_curve2
    math_piecewise_quintic_curve2
)

add_subdirectory(piecewise_curves)
