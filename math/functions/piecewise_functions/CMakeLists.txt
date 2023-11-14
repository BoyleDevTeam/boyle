add_executable(piecewise_linear_function1_test piecewise_linear_function1_test.cpp)
target_link_libraries(piecewise_linear_function1_test
  PRIVATE
    Boost::serialization
    cxxopts::cxxopts
    doctest::doctest
    Matplot++::matplot
    spdlog::spdlog
)

add_test(
  NAME piecewise_linear_function1_test
  COMMAND piecewise_linear_function1_test
)

add_executable(piecewise_cubic_function1_test piecewise_cubic_function1_test.cpp)
target_link_libraries(piecewise_cubic_function1_test
  PRIVATE
    Boost::serialization
    cxxopts::cxxopts
    doctest::doctest
    Matplot++::matplot
    spdlog::spdlog
)

add_test(
  NAME piecewise_cubic_function1_test
  COMMAND piecewise_cubic_function1_test
)
