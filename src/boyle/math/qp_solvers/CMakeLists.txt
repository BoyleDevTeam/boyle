add_header_only_library(math_qp_problem qp_problem.hpp)
target_link_libraries(math_qp_problem
  INTERFACE
    Boost::serialization
    common::logging
    math::dok_matrix
    math::lil_matrix
)
add_library(math::qp_problem ALIAS math_qp_problem)

install_libraries(
  math_qp_problem
)
