add_header_only_library(common_macros macros.hpp)
add_library(common::macros ALIAS common_macros)

add_header_only_library(common_logging logging.hpp)
target_link_libraries(common_logging
  INTERFACE
    spdlog::spdlog
)
add_library(common::logging ALIAS common_logging)

install_libraries(
  common_macros
  common_logging
)
