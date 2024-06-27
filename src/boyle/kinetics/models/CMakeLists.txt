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
    common::logging
  PUBLIC
    cvxopm::qp_problem
    cvxopm::osqp_solver
    kinetics::motion1
    kinetics::fence1
)
add_library(kinetics::route_line_cubic_acc_model ALIAS kinetics_route_line_cubic_acc_model)

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
    common::logging
  PUBLIC
    cvxopm::qp_problem
    cvxopm::osqp_solver
    kinetics::motion1
    kinetics::fence1
)
add_library(kinetics::route_line_quintic_acc_model ALIAS kinetics_route_line_quintic_acc_model)

install_libraries(
  kinetics_route_line_cubic_acc_model
  kinetics_route_line_quintic_acc_model
)
