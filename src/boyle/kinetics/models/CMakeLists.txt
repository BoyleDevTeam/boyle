add_library(kinetics_route_line_cubic_acc_model STATIC route_line_cubic_acc_model.cpp)
target_sources(kinetics_route_line_cubic_acc_model
  PUBLIC
    FILE_SET HEADERS
    FILES
      route_line_cubic_acc_model.hpp
  PRIVATE
    route_line_cubic_acc_model.cpp
)
target_link_libraries(kinetics_route_line_cubic_acc_model
  PRIVATE
    common_logging
  PUBLIC
    cvxopm_qp_problem
    cvxopm_osqp_solver
    kinetics_motion1
    kinetics_fence1
)

add_library(kinetics_route_line_quintic_acc_model STATIC route_line_quintic_acc_model.cpp)
target_sources(kinetics_route_line_quintic_acc_model
  PUBLIC
    FILE_SET HEADERS
    FILES
      route_line_quintic_acc_model.hpp
  PRIVATE
    route_line_quintic_acc_model.cpp
)
target_link_libraries(kinetics_route_line_quintic_acc_model
  PRIVATE
    common_logging
  PUBLIC
    cvxopm_qp_problem
    cvxopm_osqp_solver
    kinetics_motion1
    kinetics_fence1
)

add_library(kinetics_route_line_cubic_offset_model STATIC route_line_cubic_offset_model.cpp)
target_sources(kinetics_route_line_cubic_offset_model
  PUBLIC
    FILE_SET HEADERS
    FILES
      route_line_cubic_offset_model.hpp
  PRIVATE
    route_line_cubic_offset_model.cpp
)
target_link_libraries(kinetics_route_line_cubic_offset_model
  PRIVATE
    common_logging
  PUBLIC
    cvxopm_qp_problem
    cvxopm_osqp_solver
    kinetics_path2
    kinetics_border2
)

add_library(kinetics_route_line_quintic_offset_model STATIC route_line_quintic_offset_model.cpp)
target_sources(kinetics_route_line_quintic_offset_model
  PUBLIC
    FILE_SET HEADERS
    FILES
      route_line_quintic_offset_model.hpp
  PRIVATE
    route_line_quintic_offset_model.cpp
)
target_link_libraries(kinetics_route_line_quintic_offset_model
  PRIVATE
    common_logging
  PUBLIC
    cvxopm_qp_problem
    cvxopm_osqp_solver
    kinetics_path2
    kinetics_border2
)

install_libraries(
  kinetics_route_line_cubic_acc_model
  kinetics_route_line_quintic_acc_model
  kinetics_route_line_cubic_offset_model
  kinetics_route_line_quintic_offset_model
)
