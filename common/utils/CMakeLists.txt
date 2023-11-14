add_executable(array_view_test array_view_test.cpp)
target_link_libraries(array_view_test
  PRIVATE
    doctest::doctest
    Eigen3::Eigen
    spdlog::spdlog
)

add_test(
  NAME array_view_test
  COMMAND array_view_test
)
