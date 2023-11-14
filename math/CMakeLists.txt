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

add_subdirectory(sparse_matrix)
