add_subdirectory(models)

boyle_cxx_library(
  NAME
    kinetics_dualism
  HDRS
    "dualism.hpp"
)

boyle_cxx_library(
  NAME
    kinetics_motion
  HDRS
    "motion.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_function
)

boyle_cxx_library(
  NAME
    kinetics_route_line
  HDRS
    "route_line.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_curve
    math_utils
    math_vec2
)

boyle_cxx_library(
  NAME
    kinetics_path
  HDRS
    "path.hpp"
  DEPS
    Boost::serialization
    math_piecewise_quintic_curve
    math_utils
    math_vec2
)

boyle_cxx_library(
  NAME
    kinetics_fence
  HDRS
    "fence.hpp"
  DEPS
    Boost::serialization
    kinetics_dualism
    math_concepts
    math_piecewise_linear_function
)

boyle_cxx_library(
  NAME
    kinetics_border
  HDRS
    "border.hpp"
  DEPS
    Boost::serialization
    kinetics_dualism
    math_concepts
    math_piecewise_linear_curve
    math_vec2
)
