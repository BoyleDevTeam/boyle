add_header_only_library(cvxopm_osqp_solver osqp_solver.hpp)
target_link_libraries(cvxopm_osqp_solver
  INTERFACE
    Boost::serialization
    osqpstatic
    common_exec_on_exit
    math_csc_matrix
    cvxopm_qp_problem
)

install_libraries(
  cvxopm_osqp_solver
)
