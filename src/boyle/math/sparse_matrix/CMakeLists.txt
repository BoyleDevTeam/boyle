add_header_only_library(math_index_pair index_pair.hpp)
add_library(math::index_pair ALIAS math_index_pair)

add_header_only_library(math_dok_matrix dok_matrix.hpp)
target_link_libraries(math_dok_matrix
  INTERFACE
    Boost::serialization
    common::logging
    math::concepts
    math::index_pair
)
add_library(math::dok_matrix ALIAS math_dok_matrix)

install_libraries(
  math_index_pair
  math_dok_matrix
)
