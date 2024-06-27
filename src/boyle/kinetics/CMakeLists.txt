add_header_only_library(kinetics_dualism dualism.hpp)
add_library(kinetics::dualism ALIAS kinetics_dualism)

add_header_only_library(kinetics_motion1 motion1.hpp)
target_link_libraries(kinetics_motion1
  INTERFACE
    Boost::serialization
    math::piecewise_quintic_function1
)
add_library(kinetics::motion1 ALIAS kinetics_motion1)

add_header_only_library(kinetics_route_line2 route_line2.hpp)
target_link_libraries(kinetics_route_line2
  INTERFACE
    Boost::serialization
    math::piecewise_quintic_curve2
    math::utils
    math::vec2
)
add_library(kinetics::route_line2 ALIAS kinetics_route_line2)

add_header_only_library(kinetics_path2 path2.hpp)
target_link_libraries(kinetics_path2
  INTERFACE
    Boost::serialization
    math::piecewise_quintic_curve2
    math::utils
    math::vec2
)
add_library(kinetics::path2 ALIAS kinetics_path2)

add_header_only_library(kinetics_fence1 fence1.hpp)
target_link_libraries(kinetics_fence1
  INTERFACE
    Boost::serialization
    kinetics::dualism
    math::piecewise_linear_function1
)
add_library(kinetics::fence1 ALIAS kinetics_fence1)

add_header_only_library(kinetics_border2 border2.hpp)
target_link_libraries(kinetics_border2
  INTERFACE
    Boost::serialization
    kinetics::dualism
    math::piecewise_linear_curve2
    math::vec2
)
add_library(kinetics::border2 ALIAS kinetics_border2)

install_libraries(
  kinetics_dualism
  kinetics_motion1
  kinetics_route_line2
  kinetics_path2
  kinetics_fence1
  kinetics_border2
)
