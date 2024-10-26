add_subdirectory(models)

boyle_cxx_library(
  NAME
    kinetics_dualism
  HDRS
    "dualism.hpp"
)

boyle_cxx_library(
  NAME
    kinetics_motion1
  HDRS
    "motion1.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_function1
)

boyle_cxx_library(
  NAME
    kinetics_route_line2
  HDRS
    "route_line2.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_curve
    math_utils
    math_vec2
)

boyle_cxx_library(
  NAME
    kinetics_path2
  HDRS
    "path2.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_curve
    math_utils
    math_vec2
)

boyle_cxx_library(
  NAME
    kinetics_fence1
  HDRS
    "fence1.hpp"
  DEPS
    Boost::serialization
    kinetics_dualism
    math_concepts
    math_piecewise_linear_function1
)

boyle_cxx_library(
  NAME
    kinetics_border2
  HDRS
    "border2.hpp"
  DEPS
    Boost::serialization
    kinetics_dualism
    math_concepts
    math_piecewise_linear_curve
    math_vec2
)
