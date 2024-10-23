add_header_only_library(cvxopm_qp_problem qp_problem.hpp)
target_link_libraries(cvxopm_qp_problem
  INTERFACE
    Boost::serialization
    common_logging
    math_dok_matrix
    math_lil_matrix
)

install_libraries(
  cvxopm_qp_problem
)
