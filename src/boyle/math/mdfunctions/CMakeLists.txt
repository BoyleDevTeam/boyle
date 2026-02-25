boyle_cxx_library(
  NAME
    math_mdfunction_proxy
  HDRS
    "mdfunction_proxy.hpp"
  DEPS
    msft_proxy
)

boyle_cxx_library(
  NAME
    math_linear_mdfunction
  HDRS
    "linear_mdfunction.hpp"
  DEPS
    Boost::serialization
    math_concepts
)

boyle_cxx_library(
  NAME
    math_quadratic_mdfunction
  HDRS
    "quadratic_mdfunction.hpp"
  DEPS
    Boost::serialization
    math_concepts
)

boyle_cxx_library(
  NAME
    math_rosenbrock_function
  HDRS
    "rosenbrock_function.hpp"
  DEPS
    Boost::serialization
    math_concepts
)
