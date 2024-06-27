add_header_only_library(math_curve2_proxy curve2_proxy.hpp)
target_link_libraries(math_curve2_proxy
  INTERFACE
    msft_proxy
    math::concepts
    math::vec2
    math::utils
)
add_library(math::curve2_proxy ALIAS math_curve2_proxy)

install_libraries(
  math_curve2_proxy
)

add_subdirectory(piecewise_curves)
