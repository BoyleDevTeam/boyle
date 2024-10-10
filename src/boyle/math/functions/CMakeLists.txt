add_header_only_library(math_function1_proxy function1_proxy.hpp)
target_link_libraries(math_function1_proxy
  INTERFACE
    msft_proxy
    math_concepts
)

install_libraries(
  math_function1_proxy
)

add_subdirectory(piecewise_functions)
