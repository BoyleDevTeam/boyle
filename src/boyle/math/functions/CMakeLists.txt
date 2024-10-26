add_subdirectory(piecewise_functions)

boyle_cxx_library(
  NAME
    math_function1_proxy
  HDRS
    "function1_proxy.hpp"
  DEPS
    msft_proxy
    math_concepts
)
