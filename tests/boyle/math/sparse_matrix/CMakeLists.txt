add_boyle_test(math_dok_matrix dok_matrix_test.cpp)
target_link_libraries(math_dok_matrix_test
  PRIVATE
    math::lil_matrix
    math::coo_matrix
    math::csc_matrix
    math::csr_matrix
)

add_boyle_test(math_lil_matrix lil_matrix_test.cpp)
target_link_libraries(math_lil_matrix_test
  PRIVATE
    math::dok_matrix
    math::coo_matrix
    math::csc_matrix
    math::csr_matrix
)

add_boyle_test(math_coo_matrix coo_matrix_test.cpp)
target_link_libraries(math_coo_matrix_test
  PRIVATE
    math::dok_matrix
    math::lil_matrix
    math::csc_matrix
    math::csr_matrix
)

add_boyle_test(math_csc_matrix csc_matrix_test.cpp)
target_link_libraries(math_csc_matrix_test
  PRIVATE
    math::dok_matrix
    math::lil_matrix
    math::coo_matrix
    math::csr_matrix
)

add_boyle_test(math_csr_matrix csr_matrix_test.cpp)
target_link_libraries(math_csr_matrix_test
  PRIVATE
    math::dok_matrix
    math::lil_matrix
    math::coo_matrix
    math::csc_matrix
)
