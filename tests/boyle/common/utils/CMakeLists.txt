add_doctest_test(common_array_view array_view_test.cpp)

add_doctest_test(common_fsm fsm_test.cpp)
target_link_libraries(common_fsm_test
  PRIVATE
    common::logging
)
