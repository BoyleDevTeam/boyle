add_header_only_library(kinetics_dualism dualism.hpp)

add_header_only_library(kinetics_motion1 motion1.hpp)
target_link_libraries(kinetics_motion1
  INTERFACE
    Boost::serialization
    math_piecewise_quintic_function1
)

add_header_only_library(kinetics_route_line2 route_line2.hpp)
target_link_libraries(kinetics_route_line2
  INTERFACE
    Boost::serialization
    math_piecewise_quintic_curve2
    math_utils
    math_vec2
)

add_header_only_library(kinetics_path2 path2.hpp)
target_link_libraries(kinetics_path2
  INTERFACE
    Boost::serialization
    math_piecewise_quintic_curve2
    math_utils
    math_vec2
)

add_header_only_library(kinetics_fence1 fence1.hpp)
target_link_libraries(kinetics_fence1
  INTERFACE
    Boost::serialization
    kinetics_dualism
    math_piecewise_linear_function1
)

add_header_only_library(kinetics_border2 border2.hpp)
target_link_libraries(kinetics_border2
  INTERFACE
    Boost::serialization
    kinetics_dualism
    math_piecewise_linear_curve2
    math_vec2
)

install_libraries(
  kinetics_dualism
  kinetics_motion1
  kinetics_route_line2
  kinetics_path2
  kinetics_fence1
  kinetics_border2
)

add_subdirectory(models)
