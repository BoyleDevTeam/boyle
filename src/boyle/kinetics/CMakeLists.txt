add_header_only_library(kinetics_dualism dualism.hpp)
add_library(kinetics::dualism ALIAS kinetics_dualism)

add_header_only_library(kinetics_motion1 motion1.hpp)
target_link_libraries(kinetics_motion1
  INTERFACE
    Boost::serialization
    math::piecewise_quintic_function1
)
add_library(kinetics::motion1 ALIAS kinetics_motion1)

install_libraries(
  kinetics_dualism
  kinetics_motion1
)
