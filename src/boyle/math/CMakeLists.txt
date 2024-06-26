add_header_only_library(math_type_traits type_traits.hpp)
add_library(math::type_traits ALIAS math_type_traits)

add_header_only_library(math_concepts concepts.hpp)
add_library(math::concepts ALIAS math_concepts)

add_header_only_library(math_eigen_serialization eigen_serialization.hpp)
target_link_libraries(math_eigen_serialization
  INTERFACE
    Boost::serialization
    Eigen3::Eigen
)
add_library(math::eigen_serialization ALIAS math_eigen_serialization)

add_header_only_library(math_chebyshev chebyshev.hpp)
target_link_libraries(math_chebyshev
  INTERFACE
    Boost::serialization
    common::logging
)
add_library(math::chebyshev ALIAS math_chebyshev)

install_libraries(
  math_type_traits
  math_concepts
  math_eigen_serialization
  math_chebyshev
)
