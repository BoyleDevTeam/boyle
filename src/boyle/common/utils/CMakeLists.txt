add_header_only_library(common_macros macros.hpp)

add_header_only_library(common_logging logging.hpp)
target_link_libraries(common_logging
  INTERFACE
    spdlog::spdlog
)

add_header_only_library(common_exec_on_exit exec_on_exit.hpp)

add_header_only_library(common_fsm fsm.hpp)

install_libraries(
  common_macros
  common_logging
  common_exec_on_exit
  common_fsm
)
