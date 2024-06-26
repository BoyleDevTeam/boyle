add_header_only_library(common_macros macros.hpp)
add_library(common::macros ALIAS common_macros)

install_libraries(
  common_macros
)
