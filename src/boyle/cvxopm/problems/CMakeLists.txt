add_header_only_library(cvxopm_qp_problem qp_problem.hpp)
target_link_libraries(cvxopm_qp_problem
  INTERFACE
    Boost::serialization
    common::logging
    math::dok_matrix
    math::lil_matrix
)
add_library(cvxopm::qp_problem ALIAS cvxopm_qp_problem)

install_libraries(
  cvxopm_qp_problem
)
