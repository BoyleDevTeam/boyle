boyle_cxx_test(
  NAME
    kinetics_route_line_cubic_acc_model_test
  SRCS
    "route_line_cubic_acc_model_test.cpp"
  DEPS
    kinetics_route_line_cubic_acc_model
)

boyle_cxx_test(
  NAME
    kinetics_route_line_quintic_acc_model_test
  SRCS
    "route_line_quintic_acc_model_test.cpp"
  DEPS
    kinetics_route_line_quintic_acc_model
)

boyle_cxx_test(
  NAME
    kinetics_route_line_cubic_offset_model_test
  SRCS
    "route_line_cubic_offset_model_test.cpp"
  DEPS
    kinetics_route_line_cubic_offset_model
)

boyle_cxx_test(
  NAME
    kinetics_route_line_quintic_offset_model_test
  SRCS
    "route_line_quintic_offset_model_test.cpp"
  DEPS
    kinetics_route_line_quintic_offset_model
)
