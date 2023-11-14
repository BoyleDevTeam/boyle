add_executable(workspace_test workspace_test.cpp)
target_link_libraries(workspace_test
  PRIVATE
    doctest::doctest
)

add_test(
  NAME workspace_test
  COMMAND workspace_test
)
