add_boyle_test(math_dok_matrix dok_matrix_test.cpp)
target_link_libraries(math_dok_matrix_test
  PRIVATE
    math::lil_matrix
    math::coo_matrix
    math::csc_matrix
    math::csr_matrix
)