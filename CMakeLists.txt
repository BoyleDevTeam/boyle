cmake_minimum_required(VERSION 3.29.3)

project(boyle
  LANGUAGES
    C
    CXX
)

option(BOYLE_CHECK_PARAMS "Enable parameters checking" OFF)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_FLAGS_DEBUG "-Og -g")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")

set(CMAKE_COMPILE_WARNING_AS_ERROR True)

add_compile_options(
    -pipe
    -fno-plt
    -fexceptions
    -Wp,-D_FORTIFY_SOURCE=3
    -fstack-clash-protection
    -fcf-protection
    -fno-omit-frame-pointer
    -mno-omit-leaf-frame-pointer
    -Wall
    -Wextra
    -Wpedantic
)

if(BOYLE_CHECK_PARAMS)
    add_compile_definitions(BOYLE_CHECK_PARAMS=1)
else()
    add_compile_definitions(BOYLE_CHECK_PARAMS=0)
endif()

set(CPM_SOURCE_CACHE third_party)
set(CPM_USE_LOCAL_PACKAGES True)

include(cmake/CPM.cmake)
include(cmake/third_party.cmake)
include(cmake/utils.cmake)

include(CTest)
enable_testing()

include_directories(${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(tests)

file(CREATE_LINK
  "${CMAKE_BINARY_DIR}/compile_commands.json"
  "${CMAKE_SOURCE_DIR}/compile_commands.json"
  SYMBOLIC
)
