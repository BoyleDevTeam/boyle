add_header_only_library(math_index_pair index_pair.hpp)

add_header_only_library(math_dok_matrix dok_matrix.hpp)
target_link_libraries(math_dok_matrix
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_index_pair
)

add_header_only_library(math_lil_matrix lil_matrix.hpp)
target_link_libraries(math_lil_matrix
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_index_pair
    math_dok_matrix
)

add_header_only_library(math_coo_matrix coo_matrix.hpp)
target_link_libraries(math_coo_matrix
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_index_pair
    math_dok_matrix
)

add_header_only_library(math_csc_matrix csc_matrix.hpp)
target_link_libraries(math_csc_matrix
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_index_pair
    math_dok_matrix
)

add_header_only_library(math_csr_matrix csr_matrix.hpp)
target_link_libraries(math_csr_matrix
  INTERFACE
    Boost::serialization
    common_logging
    math_concepts
    math_index_pair
    math_dok_matrix
)

add_header_only_library(math_sparse_matrix_proxy sparse_matrix_proxy.hpp)
target_link_libraries(math_sparse_matrix_proxy
  INTERFACE
    msft_proxy
    math_concepts
)

install_libraries(
  math_index_pair
  math_dok_matrix
  math_lil_matrix
  math_coo_matrix
  math_csc_matrix
  math_csr_matrix
  math_sparse_matrix_proxy
)
