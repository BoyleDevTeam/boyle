add_executable(piecewise_linear_curve2_test piecewise_linear_curve2_test.cpp)
target_link_libraries(piecewise_linear_curve2_test
  PRIVATE
    Boost::serialization
    cxxopts::cxxopts
    doctest::doctest
    Matplot++::matplot
    spdlog::spdlog
)

add_test(
  NAME piecewise_linear_curve2_test
  COMMAND piecewise_linear_curve2_test
)
