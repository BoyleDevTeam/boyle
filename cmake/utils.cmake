function(add_header_only_library name header_name)
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE LIBRARY_HEADER_ONLY)
  target_sources(${name}
    INTERFACE
      FILE_SET HEADERS
      FILES
        ${header_name}
  )
endfunction()

function(install_libraries)
  file(RELATIVE_PATH current_subdirectory "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
  install(
    TARGETS
      ${ARGV}
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${current_subdirectory}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  unset(current_subdirectory)
endfunction()

function(add_doctest_test name source_name)
  set(test_name ${name}_test)
  add_executable(${test_name} ${source_name})
  target_link_libraries(${test_name}
    PRIVATE
      cxxopts::cxxopts
      doctest::doctest
      Matplot++::matplot
      ${name}
  )
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
  )
endfunction()
