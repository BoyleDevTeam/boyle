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

install_libraries(
  math_type_traits
  math_concepts
  math_array
)
