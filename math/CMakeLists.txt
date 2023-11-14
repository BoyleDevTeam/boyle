add_executable(vec2_test vec2_test.cpp)
target_link_libraries(vec2_test
  PRIVATE
    doctest::doctest
    Boost::serialization
)

add_test(
  NAME vec2_test
  COMMAND vec2_test
)

add_executable(vec3_test vec3_test.cpp)
target_link_libraries(vec3_test
  PRIVATE
    doctest::doctest
    Boost::serialization
)

add_test(
  NAME vec3_test
  COMMAND vec3_test
)

add_executable(cubic_interpolation_test cubic_interpolation_test.cpp)
target_link_libraries(cubic_interpolation_test
  PRIVATE
    doctest::doctest
    spdlog::spdlog
)

add_test(
  NAME cubic_interpolation_test
  COMMAND cubic_interpolation_test
)

add_executable(quintic_interpolation_test quintic_interpolation_test.cpp)
target_link_libraries(quintic_interpolation_test
  PRIVATE
    doctest::doctest
    spdlog::spdlog
)

add_test(
  NAME quintic_interpolation_test
  COMMAND quintic_interpolation_test
)

add_subdirectory(functions)
add_subdirectory(qp_solvers)
add_subdirectory(sparse_matrix)
