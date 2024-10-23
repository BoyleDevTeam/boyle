cmake_minimum_required(VERSION 3.30.0)

project(Boyle
  LANGUAGES
    C
    CXX
  VERSION
    0.1.0
)

option(BOYLE_CHECK_PARAMS "Enable parameters checking" OFF)

if(BOYLE_CHECK_PARAMS)
  add_compile_definitions(BOYLE_CHECK_PARAMS=1)
else()
  add_compile_definitions(BOYLE_CHECK_PARAMS=0)
endif()

set(CPM_SOURCE_CACHE third_party)
set(CPM_USE_LOCAL_PACKAGES True)

include(cmake/utils.cmake)
include(cmake/CPM.cmake)
include(cmake/third_party.cmake)

include(CTest)
enable_testing()

include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(tests)

include(cmake/export.cmake)

if(CMAKE_GENERATOR STREQUAL "Ninja")
  file(CREATE_LINK
    "${CMAKE_BINARY_DIR}/compile_commands.json"
    "${CMAKE_SOURCE_DIR}/compile_commands.json"
    SYMBOLIC
  )
endif()
