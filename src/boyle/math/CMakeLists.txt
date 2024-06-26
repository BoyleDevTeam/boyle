add_header_only_library(math_type_traits type_traits.hpp)
add_library(math::type_traits ALIAS math_type_traits)

add_header_only_library(math_concepts concepts.hpp)
add_library(math::concepts ALIAS math_concepts)

add_header_only_library(math_array array.hpp)
target_link_libraries(math_array
  INTERFACE
    Boost::serialization
    OpenBLAS
    LPACKE
)
add_library(math::array ALIAS math_array)

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

install_libraries(
  math_type_traits
  math_concepts
  math_array
  math_chebyshev
  math_duplet
  math_vec2
  math_triplet
)
