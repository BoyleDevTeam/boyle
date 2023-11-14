add_library(route_line_builder STATIC route_line_builder.cpp)
target_sources(route_line_builder
  PUBLIC
    route_line_builder.h
    ${CMAKE_SOURCE_DIR}/tasks/task.hpp
    ${CMAKE_SOURCE_DIR}/workspaces/workspace.hpp
  INTERFACE
    ${CMAKE_SOURCE_DIR}/common/utils/macros.hpp
  PRIVATE
    route_line_builder.cpp
)

add_library(path_optimizer STATIC path_optimizer.cpp)
target_sources(path_optimizer
  PUBLIC
    path_optimizer.h
    ${CMAKE_SOURCE_DIR}/tasks/task.hpp
    ${CMAKE_SOURCE_DIR}/workspaces/workspace.hpp
  INTERFACE
    ${CMAKE_SOURCE_DIR}/common/utils/macros.hpp
  PRIVATE
    path_optimizer.cpp
)

add_library(speed_optimizer STATIC speed_optimizer.cpp)
target_sources(speed_optimizer
  PUBLIC
    speed_optimizer.h
    ${CMAKE_SOURCE_DIR}/tasks/task.hpp
    ${CMAKE_SOURCE_DIR}/workspaces/workspace.hpp
  INTERFACE
    ${CMAKE_SOURCE_DIR}/common/utils/macros.hpp
  PRIVATE
    speed_optimizer.cpp
)
