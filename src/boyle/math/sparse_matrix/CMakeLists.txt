add_header_only_library(math_index_pair index_pair.hpp)
add_library(math::index_pair ALIAS math_index_pair)

install_libraries(
  math_index_pair
)
