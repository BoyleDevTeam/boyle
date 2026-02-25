add_subdirectory(curves)
add_subdirectory(functions)
add_subdirectory(mdfunctions)
add_subdirectory(dense)
add_subdirectory(sparse)

boyle_cxx_library(
  NAME
    math_type_traits
  HDRS
    "type_traits.hpp"
)

boyle_cxx_library(
  NAME
    math_concepts
  HDRS
    "concepts.hpp"
)

boyle_cxx_library(
  NAME
    math_fft
  HDRS
    "fft.hpp"
  DEPS
    pocketfft
)

boyle_cxx_library(
  NAME
    math_chebyshev
  HDRS
    "chebyshev.hpp"
  DEPS
    Boost::serialization
)

boyle_cxx_library(
  NAME
    math_utils
  HDRS
    "utils.hpp"
  DEPS
    math_concepts
)

boyle_cxx_library(
  NAME
    math_cubic_interpolation
  HDRS
    "cubic_interpolation.hpp"
  DEPS
    math_concepts
    math_utils
)

boyle_cxx_library(
  NAME
    math_quintic_interpolation
  HDRS
    "quintic_interpolation.hpp"
  DEPS
    math_concepts
    math_cubic_interpolation
)
