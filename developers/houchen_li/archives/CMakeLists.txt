boyle_cxx_library(
  NAME
    math_eigen_serialization
  HDRS
    "eigen_serialization.hpp"
  DEPS
    Boost::serialization
    Eigen3::Eigen
)

boyle_cxx_test(
  NAME
    math_eigen_serialization_test
  SRCS
    "eigen_serialization_test.cpp"
  DEPS
    math_eigen_serialization
)

boyle_cxx_library(
  NAME
    math_type_concepts
  HDRS
    "type_concepts.hpp"
  DEPS
    math_type_traits
)
