add_header_only_library(math_eigen_serialization eigen_serialization.hpp)
target_link_libraries(math_eigen_serialization
  INTERFACE
    Boost::serialization
    Eigen3::Eigen
)

install_libraries(
  math_eigen_serialization
)
