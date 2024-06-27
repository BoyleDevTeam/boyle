add_header_only_library(cvxopm_osqp_solver osqp_solver.hpp)
target_link_libraries(cvxopm_osqp_solver
  INTERFACE
    Boost::serialization
    osqpstatic
    common::exec_on_exit
    math::csc_matrix
    cvxopm::qp_problem
)
add_library(cvxopm::osqp_solver ALIAS cvxopm_osqp_solver)

install_libraries(
  cvxopm_osqp_solver
)
