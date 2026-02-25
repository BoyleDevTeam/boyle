boyle_cxx_library(
  NAME
    math_dense_traits
  HDRS
    "dense_traits.hpp"
  DEPS
    common_aligned_allocator
    math_concepts
)

boyle_cxx_library(
  NAME
    math_vector_view
  HDRS
    "vector_view.hpp"
  DEPS
    math_concepts
    math_dense_traits
)

boyle_cxx_library(
  NAME
    math_vector
  HDRS
    "vector.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
    math_dense_traits
    math_vector_view
)

boyle_cxx_library(
  NAME
    math_vectorx
  HDRS
    "vectorx.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
    math_dense_traits
    math_vector_view
)

boyle_cxx_library(
  NAME
    math_matrix_view
  HDRS
    "matrix_view.hpp"
  DEPS
    math_concepts
    math_dense_traits
)

boyle_cxx_library(
  NAME
    math_matrix
  HDRS
    "matrix.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
    math_dense_traits
    math_matrix_view
)

boyle_cxx_library(
  NAME
    math_matrixx
  HDRS
    "matrixx.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
    math_dense_traits
    math_matrix_view
)

boyle_cxx_library(
  NAME
    math_lu_dcmp
  HDRS
    "lu_dcmp.hpp"
    "detail/dense_partial_pivot_trait.hpp"
  DEPS
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
)

boyle_cxx_library(
  NAME
    math_qr_dcmp
  HDRS
    "qr_dcmp.hpp"
    "detail/dense_degenerate_trait.hpp"
    "detail/dense_generate_trait.hpp"
  DEPS
    "$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,blas_lapack,blas_lapack_static>"
    math_concepts
)

boyle_cxx_library(
  NAME
    math_vec2
  HDRS
    "vec2.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    math_concepts
    math_dense_traits
)

boyle_cxx_library(
  NAME
    math_vec3
  HDRS
    "vec3.hpp"
    "detail/dense_norm_trait.hpp"
  DEPS
    Boost::serialization
    common_in_in_in_out_result
    math_concepts
    math_dense_traits
)

if(BOYLE_USE_SIMD STREQUAL "AVX512")
  set_target_properties(math_vec2 math_vec3 PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "BOYLE_USE_SIMD=1"
    INTERFACE_COMPILE_OPTIONS "-mavx512f"
  )
endif()
