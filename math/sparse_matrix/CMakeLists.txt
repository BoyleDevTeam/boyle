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
