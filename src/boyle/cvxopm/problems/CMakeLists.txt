boyle_cxx_library(
  NAME
    cvxopm_qp_problem
  HDRS
    "qp_problem.hpp"
  DEPS
    Boost::serialization
    common_aligned_allocator
    math_dok_vector
    math_dok_matrix
    math_lil_matrix
)
