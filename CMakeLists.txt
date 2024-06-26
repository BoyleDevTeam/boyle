cmake_minimum_required(VERSION 3.30.0)

project(boyle
  LANGUAGES
    C
    CXX
)

option(BOYLE_CHECK_PARAMS "Enable parameters checking" OFF)

if(BOYLE_CHECK_PARAMS)
    add_compile_definitions(BOYLE_CHECK_PARAMS=1)
else()
    add_compile_definitions(BOYLE_CHECK_PARAMS=0)
endif()

set(CPM_SOURCE_CACHE third_party)
set(CPM_USE_LOCAL_PACKAGES True)

include(cmake/CPM.cmake)
include(cmake/utils.cmake)

include(CTest)
enable_testing()

file(CREATE_LINK
  "${CMAKE_BINARY_DIR}/compile_commands.json"
  "${CMAKE_SOURCE_DIR}/compile_commands.json"
  SYMBOLIC
)
