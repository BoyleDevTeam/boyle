add_executable(qp_problem_test qp_problem_test.cpp)
target_link_libraries(qp_problem_test
  PRIVATE
    Boost::serialization
    doctest::doctest
    osqpstatic
    spdlog::spdlog
)

add_test(
  NAME qp_problem_test
  COMMAND qp_problem_test
)
