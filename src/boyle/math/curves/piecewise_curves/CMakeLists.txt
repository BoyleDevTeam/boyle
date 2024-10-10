add_header_only_library(math_piecewise_linear_curve2 piecewise_linear_curve2.hpp)
target_link_libraries(math_piecewise_linear_curve2
  INTERFACE 
    Boost::serialization
    common_logging
    math_duplet
    math_piecewise_linear_function1
    math_utils
    math_vec2
)

add_header_only_library(math_piecewise_cubic_curve2 piecewise_cubic_curve2.hpp)
target_link_libraries(math_piecewise_cubic_curve2
  INTERFACE
    Boost::serialization
    common_logging
    math_duplet
    math_cubic_interpolation
    math_piecewise_cubic_function1
    math_utils
    math_vec2
)

add_header_only_library(math_piecewise_quintic_curve2 piecewise_quintic_curve2.hpp)
target_link_libraries(math_piecewise_quintic_curve2
  INTERFACE
    Boost::serialization
    common_logging
    math_duplet
    math_piecewise_quintic_function1
    math_quintic_interpolation
    math_utils
    math_vec2
)

install_libraries(
  math_piecewise_linear_curve2
  math_piecewise_cubic_curve2
  math_piecewise_quintic_curve2
)
