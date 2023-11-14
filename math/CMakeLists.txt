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
