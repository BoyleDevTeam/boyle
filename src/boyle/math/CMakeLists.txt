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

add_header_only_library(math_duplet duplet.hpp)
target_link_libraries(math_duplet
  INTERFACE
    Boost::serialization
)
add_library(math::duplet ALIAS math_duplet)

add_header_only_library(math_vec2 vec2.hpp)
target_link_libraries(math_vec2
  INTERFACE
    math::concepts
    math::duplet
)
add_library(math::vec2 ALIAS math_vec2)

add_header_only_library(math_triplet triplet.hpp)
target_link_libraries(math_triplet
  INTERFACE
    Boost::serialization
)
add_library(math::triplet ALIAS math_triplet)

add_header_only_library(math_vec3 vec3.hpp)
target_link_libraries(math_vec3
  INTERFACE
    math::concepts
    math::triplet
)
add_library(math::vec3 ALIAS math_vec3)

add_header_only_library(math_utils utils.hpp)
target_link_libraries(math_utils
  INTERFACE
    common::logging
    math::concepts
    math::duplet
    math::triplet
    math::vec2
    math::vec3
)
add_library(math::utils ALIAS math_utils)

add_header_only_library(math_cubic_interpolation cubic_interpolation.hpp)
target_link_libraries(math_cubic_interpolation
  INTERFACE
    math::concepts
    math::utils
)
add_library(math::cubic_interpolation ALIAS math_cubic_interpolation)

add_header_only_library(math_quintic_interpolation quintic_interpolation.hpp)
target_link_libraries(math_quintic_interpolation
  INTERFACE
    math::concepts
    math::cubic_interpolation
)
add_library(math::quintic_interpolation ALIAS math_quintic_interpolation)

install_libraries(
  math_type_traits
  math_concepts
  math_eigen_serialization
  math_chebyshev
  math_duplet
  math_vec2
  math_triplet
  math_vec3
  math_utils
  math_cubic_interpolation
  math_quintic_interpolation
)

add_subdirectory(curves)
add_subdirectory(functions)
add_subdirectory(sparse_matrix)
