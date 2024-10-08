add_header_only_library(math_type_traits type_traits.hpp)

add_header_only_library(math_concepts concepts.hpp)

add_header_only_library(math_array array.hpp)
target_link_libraries(math_array
  INTERFACE
    Boost::serialization
    LPACKE
)

add_header_only_library(math_chebyshev chebyshev.hpp)
target_link_libraries(math_chebyshev
  INTERFACE
    Boost::serialization
    common_logging
)

add_header_only_library(math_duplet duplet.hpp)
target_link_libraries(math_duplet
  INTERFACE
    Boost::serialization
)

add_header_only_library(math_vec2 vec2.hpp)
target_link_libraries(math_vec2
  INTERFACE
    math_concepts
    math_duplet
)

add_header_only_library(math_triplet triplet.hpp)
target_link_libraries(math_triplet
  INTERFACE
    Boost::serialization
)

add_header_only_library(math_vec3 vec3.hpp)
target_link_libraries(math_vec3
  INTERFACE
    math_concepts
    math_triplet
)

add_header_only_library(math_utils utils.hpp)
target_link_libraries(math_utils
  INTERFACE
    common_logging
    math_concepts
    math_duplet
    math_triplet
    math_vec2
    math_vec3
)

add_header_only_library(math_cubic_interpolation cubic_interpolation.hpp)
target_link_libraries(math_cubic_interpolation
  INTERFACE
    math_concepts
    math_utils
)

add_header_only_library(math_quintic_interpolation quintic_interpolation.hpp)
target_link_libraries(math_quintic_interpolation
  INTERFACE
    math_concepts
    math_cubic_interpolation
)

install_libraries(
  math_type_traits
  math_concepts
  math_array
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
