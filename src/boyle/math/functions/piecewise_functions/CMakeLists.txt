add_header_only_library(math_piecewise_linear_function1 piecewise_linear_function1.hpp)
target_link_libraries(math_piecewise_linear_function1
  INTERFACE
    Boost::serialization
    common::logging
    math::concepts
    math::utils
)
add_library(math::piecewise_linear_function1 ALIAS math_piecewise_linear_function1)

add_header_only_library(math_piecewise_cubic_function1 piecewise_cubic_function1.hpp)
target_link_libraries(math_piecewise_cubic_function1
  INTERFACE
    Boost::serialization
    common::logging
    math::concepts
    math::cubic_interpolation
)
add_library(math::piecewise_cubic_function1 ALIAS math_piecewise_cubic_function1)

add_header_only_library(math_piecewise_quintic_function1 piecewise_quintic_function1.hpp)
target_link_libraries(math_piecewise_quintic_function1
  INTERFACE
    Boost::serialization
    common::logging
    math::concepts
    math::quintic_interpolation
)
add_library(math::piecewise_quintic_function1 ALIAS math_piecewise_quintic_function1)

install_libraries(
  math_piecewise_linear_function1
  math_piecewise_cubic_function1
  math_piecewise_quintic_function1
)
