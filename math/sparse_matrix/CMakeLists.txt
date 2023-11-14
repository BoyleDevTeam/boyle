add_executable(coo_matrix_test coo_matrix_test.cpp)
target_link_libraries(coo_matrix_test
  PRIVATE
    Boost::serialization
    doctest::doctest
    spdlog::spdlog
)

add_test(
  NAME coo_matrix_test
  COMMAND coo_matrix_test
)

add_executable(lil_matrix_test lil_matrix_test.cpp)
target_link_libraries(lil_matrix_test
  PRIVATE
    Boost::serialization
    doctest::doctest
    spdlog::spdlog
)

add_test(
  NAME lil_matrix_test
  COMMAND lil_matrix_test
)
