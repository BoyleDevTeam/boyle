add_header_only_library(common_macros macros.hpp)
add_library(common::macros ALIAS common_macros)

add_header_only_library(common_logging logging.hpp)
target_link_libraries(common_logging
  INTERFACE
    spdlog::spdlog
)
add_library(common::logging ALIAS common_logging)

add_header_only_library(common_exec_on_exit exec_on_exit.hpp)
add_library(common::exec_on_exit ALIAS common_exec_on_exit)

add_header_only_library(common_array_view array_view.hpp)
target_link_libraries(common_array_view
  INTERFACE
    Eigen3::Eigen
)
add_library(common::array_view ALIAS common_array_view)

add_header_only_library(common_fsm fsm.hpp)
add_library(common::fsm ALIAS common_fsm)

install_libraries(
  common_macros
  common_logging
  common_exec_on_exit
  common_array_view
  common_fsm
)
