add_header_only_library(math_piecewise_linear_function1 piecewise_linear_function1.hpp)
target_link_libraries(math_piecewise_linear_function1
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_utils
)

add_header_only_library(math_piecewise_cubic_function1 piecewise_cubic_function1.hpp)
target_link_libraries(math_piecewise_cubic_function1
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_cubic_interpolation
)

add_header_only_library(math_piecewise_quintic_function1 piecewise_quintic_function1.hpp)
target_link_libraries(math_piecewise_quintic_function1
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_quintic_interpolation
)

install_libraries(
  math_piecewise_linear_function1
  math_piecewise_cubic_function1
  math_piecewise_quintic_function1
)
