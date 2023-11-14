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

add_executable(osqp_solver_test osqp_solver_test.cpp)
target_link_libraries(osqp_solver_test
  PRIVATE
    Boost::serialization
    doctest::doctest
    osqpstatic
    spdlog::spdlog
)

add_test(
  NAME osqp_solver_test
  COMMAND osqp_solver_test
)
