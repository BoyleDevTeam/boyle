add_subdirectory(problems)
add_subdirectory(solvers)

boyle_cxx_library(
  NAME
    cvxopm_settings
  HDRS
    "settings.hpp"
  DEPS
    Boost::serialization
    common_aligned_memory_resource
)

boyle_cxx_library(
  NAME
    cvxopm_result
  HDRS
    "result.hpp"
  DEPS
    Boost::serialization
)

boyle_cxx_library(
  NAME
    cvxopm_info
  HDRS
    "info.hpp"
  DEPS
    Boost::serialization
)
