boyle_cxx_library(
  NAME
    math_function1_proxy
  HDRS
    "function1_proxy.hpp"
  DEPS
    msft_proxy
    math_concepts
)

boyle_cxx_library(
  NAME
    math_piecewise_linear_function1
  HDRS
    "piecewise_linear_function1.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_utils
    math_vec2
    math_vec3
)

boyle_cxx_library(
  NAME
    math_piecewise_cubic_function1
  HDRS
    "piecewise_cubic_function1.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_cubic_interpolation
    math_vec2
    math_vec3
)

boyle_cxx_library(
  NAME
    math_piecewise_quintic_function1
  HDRS
    "piecewise_quintic_function1.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_quintic_interpolation
    math_vec2
    math_vec3
)
