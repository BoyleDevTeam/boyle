add_library(route_line_acc_model STATIC route_line_acc_model.cpp)
target_sources(route_line_acc_model
  PUBLIC
    route_line_acc_model.h
  PRIVATE
    route_line_acc_model.cpp
)
target_link_libraries(route_line_acc_model
  PUBLIC
    Boost::serialization
    osqpstatic
    spdlog::spdlog
)

add_executable(route_line_acc_model_test route_line_acc_model_test.cpp)
target_link_libraries(route_line_acc_model_test
  PRIVATE
    cxxopts::cxxopts
    doctest::doctest
    Matplot++::matplot
    route_line_acc_model
    spdlog::spdlog
)

add_test(
  NAME route_line_acc_model_test
  COMMAND route_line_acc_model_test
)
