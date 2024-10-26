add_custom_command(
  OUTPUT "${CMAKE_SOURCE_DIR}/compile_commands.json"
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          "${CMAKE_BINARY_DIR}/compile_commands.json"
          "${CMAKE_SOURCE_DIR}/compile_commands.json"
  DEPENDS "${CMAKE_BINARY_DIR}/compile_commands.json"
  COMMENT "Create symlink for compile_commands.json to project root"
  VERBATIM
)

add_custom_target(compile_commands ALL DEPENDS "${CMAKE_SOURCE_DIR}/compile_commands.json")
