add_header_only_library(kinetics_dualism dualism.hpp)
add_library(kinetics::dualism ALIAS kinetics_dualism)

install_libraries(
  kinetics_dualism
)
