add_header_only_library(math_curve2_proxy curve2_proxy.hpp)
target_link_libraries(math_curve2_proxy
  INTERFACE
    msft_proxy
    math_concepts
    math_vec2
    math_utils
)

install_libraries(
  math_curve2_proxy
)

add_subdirectory(piecewise_curves)
