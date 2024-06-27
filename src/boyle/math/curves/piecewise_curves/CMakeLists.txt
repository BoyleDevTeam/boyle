add_header_only_library(math_piecewise_linear_curve2 piecewise_linear_curve2.hpp)
target_link_libraries(math_piecewise_linear_curve2
  INTERFACE 
    Boost::serialization
    common::logging
    math::duplet
    math::piecewise_linear_function1
    math::utils
    math::vec2
)
add_library(math::piecewise_linear_curve2 ALIAS math_piecewise_linear_curve2)

add_header_only_library(math_piecewise_cubic_curve2 piecewise_cubic_curve2.hpp)
target_link_libraries(math_piecewise_cubic_curve2
  INTERFACE
    Boost::serialization
    common::logging
    math::duplet
    math::cubic_interpolation
    math::piecewise_cubic_function1
    math::utils
    math::vec2
)
add_library(math::piecewise_cubic_curve2 ALIAS math_piecewise_cubic_curve2)

add_header_only_library(math_piecewise_quintic_curve2 piecewise_quintic_curve2.hpp)
target_link_libraries(math_piecewise_quintic_curve2
  INTERFACE
    Boost::serialization
    common::logging
    math::duplet
    math::piecewise_quintic_function1
    math::quintic_interpolation
    math::utils
    math::vec2
)
add_library(math::piecewise_quintic_curve2 ALIAS math_piecewise_quintic_curve2)

install_libraries(
  math_piecewise_linear_curve2
  math_piecewise_cubic_curve2
  math_piecewise_quintic_curve2
)
