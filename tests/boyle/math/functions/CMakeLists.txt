add_boyle_test(math_function1_proxy function1_proxy_test.cpp)
target_link_libraries(math_function1_proxy_test
  PRIVATE
    math::piecewise_linear_function1
    math::piecewise_cubic_function1
    math::piecewise_quintic_function1
)

add_subdirectory(piecewise_functions)
