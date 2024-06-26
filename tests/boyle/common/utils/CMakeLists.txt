add_boyle_test(common_fsm fsm_test.cpp)
target_link_libraries(common_fsm_test
  PRIVATE
    common::logging
)
