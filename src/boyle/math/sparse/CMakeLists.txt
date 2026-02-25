boyle_cxx_library(
  NAME
    math_index_pair
  HDRS
    "index_pair.hpp"
)

boyle_cxx_library(
  NAME
    math_sparse_traits
  HDRS
    "sparse_traits.hpp"
  DEPS
    common_aligned_allocator
    math_concepts
    math_index_pair
)

boyle_cxx_library(
  NAME
    math_sparse_matrix_proxy
  HDRS
    "sparse_matrix_proxy.hpp"
  DEPS
    msft_proxy
    math_concepts
)

boyle_cxx_library(
  NAME
    math_dok_matrix
  HDRS
    "dok_matrix.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_sparse_traits
)

boyle_cxx_library(
  NAME
    math_lil_matrix
  HDRS
    "lil_matrix.hpp"
  DEPS
    Boost::serialization
    common_aligned_allocator
    math_concepts
    math_sparse_traits
)

boyle_cxx_library(
  NAME
    math_coo_matrix
  HDRS
    "coo_matrix.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_sparse_traits
)

boyle_cxx_library(
  NAME
    math_csc_matrix
  HDRS
    "csc_matrix.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_index_pair
    math_sparse_traits
)

boyle_cxx_library(
  NAME
    math_csr_matrix
  HDRS
    "csr_matrix.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_index_pair
    math_sparse_traits
)

if(BOYLE_USE_BOOST_UNORDERED)
  set_target_properties(math_sparse_traits math_index_pair math_dok_matrix math_lil_matrix math_coo_matrix math_csc_matrix math_csr_matrix
    PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_BOOST_UNORDERED=1"
  )
endif()
