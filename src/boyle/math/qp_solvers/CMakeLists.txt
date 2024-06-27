add_header_only_library(math_qp_problem qp_problem.hpp)
target_link_libraries(math_qp_problem
  INTERFACE
    Boost::serialization
    common::logging
    math::dok_matrix
    math::lil_matrix
)
add_library(math::qp_problem ALIAS math_qp_problem)

add_header_only_library(math_osqp_solver osqp_solver.hpp)
target_link_libraries(math_osqp_solver
  INTERFACE
    Boost::serialization
    osqpstatic
    common::exec_on_exit
    math::csc_matrix
    math::qp_problem
)
add_library(math::osqp_solver ALIAS math_osqp_solver)

install_libraries(
  math_qp_problem
  math_osqp_solver
)
