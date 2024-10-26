add_subdirectory(piecewise_curves)

boyle_cxx_library(
  NAME
    math_curve2_proxy
  HDRS
    "curve2_proxy.hpp"
  DEPS
    msft_proxy
    math_concepts
    math_duplet
    math_vec2
    math_utils
)
