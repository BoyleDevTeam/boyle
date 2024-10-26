include(CMakeParseArguments)

if(BOYLE_BUILD_TESTING)
  include(CTest)
  enable_testing()
endif()

if(NOT DEFINED BOYLE_IDE_FOLDER)
  set(BOYLE_IDE_FOLDER ${CMAKE_PROJECT_NAME})
endif()

function(boyle_cxx_library)
  cmake_parse_arguments(BOYLE_CXX_LIB
    "DISABLE_INSTALL;PUBLIC;TESTONLY"
    "NAME"
    "HDRS;SRCS;COPTS;DEFINES;LINKOPTS;DEPS"
    ${ARGN}
  )

  set(_SRCS "${BOYLE_CXX_LIB_SRCS}")
  list(FILTER _SRCS EXCLUDE REGEX ".*\\.(h|hpp|inc)")
  if(_SRCS STREQUAL "")
    set(_IS_INTERFACE 1)
  else()
    set(_IS_INTERFACE 0)
  endif()
  unset(_SRCS)

  if(NOT _IS_INTERFACE)
    add_library(${BOYLE_CXX_LIB_NAME} "")
    set_property(TARGET ${BOYLE_CXX_LIB_NAME} PROPERTY LINKER_LANGUAGE "CXX")
    target_sources(${BOYLE_CXX_LIB_NAME}
      PRIVATE
        ${BOYLE_CXX_LIB_SRCS}
      PUBLIC
        FILE_SET HEADERS
          TYPE HEADERS
          BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
          FILES ${BOYLE_CXX_LIB_HDRS}
    )
    target_include_directories(${BOYLE_CXX_LIB_NAME}
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    target_compile_options(${BOYLE_CXX_LIB_NAME} PRIVATE ${BOYLE_CXX_LIB_COPTS})
    target_compile_definitions(${BOYLE_CXX_LIB_NAME} PUBLIC ${BOYLE_CXX_LIB_DEFINES})
    target_link_options(${BOYLE_CXX_LIB_NAME} PRIVATE ${BOYLE_CXX_LIB_LINKOPTS})
    target_link_libraries(${BOYLE_CXX_LIB_NAME} PUBLIC ${BOYLE_CXX_LIB_DEPS})

    if(BOYLE_CXX_LIB_PUBLIC)
      set_property(TARGET ${BOYLE_CXX_LIB_NAME} PROPERTY FOLDER ${BOYLE_IDE_FOLDER})
    elseif(BOYLE_CXX_LIB_TESTONLY)
      set_property(TARGET ${BOYLE_CXX_LIB_NAME} PROPERTY FOLDER ${BOYLE_IDE_FOLDER}/tests)
    else()
      set_property(TARGET ${BOYLE_CXX_LIB_NAME} PROPERTY FOLDER ${BOYLE_IDE_FOLDER}/internal)
    endif()
    
    if(BOYLE_ENABLE_INSTALL)
      set_target_properties(${BOYLE_CXX_LIB_NAME} PROPERTIES
        OUTPUT_NAME "boyle_${BOYLE_CXX_LIB_NAME}"
        SOVERSION "${BOYLE_SOVERSION}"
      )
    endif()
  else()
    add_library(${BOYLE_CXX_LIB_NAME} INTERFACE)
    target_sources(${BOYLE_CXX_LIB_NAME}
      INTERFACE
        FILE_SET HEADERS
          TYPE HEADERS
          BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
          FILES ${BOYLE_CXX_LIB_HDRS}
    )
    target_include_directories(${BOYLE_CXX_LIB_NAME}
      INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    target_compile_definitions(${BOYLE_CXX_LIB_NAME} INTERFACE ${BOYLE_CXX_LIB_DEFINES})
    target_link_options(${BOYLE_CXX_LIB_NAME} INTERFACE ${BOYLE_CXX_LIB_LINKOPTS})
    target_link_libraries(${BOYLE_CXX_LIB_NAME} INTERFACE ${BOYLE_CXX_LIB_DEPS})
  endif()

  unset(_IS_INTERFACE)

  if(BOYLE_ENABLE_INSTALL)
    file(RELATIVE_PATH _REL_PATH "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
    install(TARGETS ${BOYLE_CXX_LIB_NAME} EXPORT ${CMAKE_PROJECT_NAME}Targets
      FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${_REL_PATH}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    unset(_REL_PATH)
  endif()

  add_library(${CMAKE_PROJECT_NAME}::${BOYLE_CXX_LIB_NAME} ALIAS ${BOYLE_CXX_LIB_NAME})
endfunction()

function(boyle_cxx_module)
  cmake_parse_arguments(BOYLE_CXX_MODULE
    "DISABLE_INSTALL;PUBLIC;TESTONLY"
    "NAME"
    "IXXS;SRCS;COPTS;DEFINES;LINKOPTS;DEPS"
    ${ARGN}
  )

  add_library(${BOYLE_CXX_MODULE_NAME} "")
  target_sources(${BOYLE_CXX_MODULE_NAME}
    PUBLIC
      FILE_SET CXX_MODULES
        TYPE CXX_MODULES
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        FILES ${BOYLE_CXX_MODULE_IXXS}
  )
  target_include_directories(${BOYLE_CXX_MODULE_NAME}
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/modules>"
  )
  target_compile_options(${BOYLE_CXX_MODULE_NAME} PRIVATE ${BOYLE_CXX_MODULE_COPTS})
  target_compile_definitions(${BOYLE_CXX_MODULE_NAME} PUBLIC ${BOYLE_CXX_MODULE_DEFINES})
  target_link_options(${BOYLE_CXX_MODULE_NAME} PRIVATE ${BOYLE_CXX_MODULE_LINKOPTS})
  target_link_libraries(${BOYLE_CXX_MODULE_NAME} PUBLIC ${BOYLE_CXX_MODULE_DEPS})

  if(BOYLE_ENABLE_INSTALL)
    file(RELATIVE_PATH _REL_PATH "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
    install(TARGETS ${BOYLE_CXX_MODULE_NAME} EXPORT ${CMAKE_PROJECT_NAME}Targets
      FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_PREFIX}/modules/${_REL_PATH}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    unset(_REL_PATH)
  endif()

  add_library(${CMAKE_PROJECT_NAME}::${BOYLE_CXX_MODULE_NAME} ALIAS ${BOYLE_CXX_MODULE_NAME})
endfunction()

function(boyle_cxx_test)
  if(NOT (BUILD_TESTING AND BOYLE_BUILD_TESTING))
    return()
  endif()

  cmake_parse_arguments(BOYLE_CXX_TEST
    ""
    "NAME"
    "SRCS;COPTS;DEFINES;LINKOPTS;DEPS"
    ${ARGN}
  )

  add_executable(${BOYLE_CXX_TEST_NAME} ${BOYLE_CXX_TEST_SRCS})
  set_property(TARGET ${BOYLE_CXX_TEST_NAME} PROPERTY FOLDER ${BOYLE_IDE_FOLDER}/tests)
  target_include_directories(${BOYLE_CXX_TEST_NAME}
    PRIVATE
      "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  target_compile_options(${BOYLE_CXX_TEST_NAME} PRIVATE ${BOYLE_CXX_TEST_COPTS})
  target_compile_definitions(${BOYLE_CXX_TEST_NAME} PRIVATE ${BOYLE_CXX_TEST_DEFINES})
  target_link_options(${BOYLE_CXX_TEST_NAME} PRIVATE ${BOYLE_CXX_TEST_LINKOPTS})
  target_link_libraries(${BOYLE_CXX_TEST_NAME}
    PRIVATE
      cxxopts::cxxopts
      doctest::doctest
      Matplot++::matplot
      ${BOYLE_CXX_TEST_DEPS}
  )

  add_test(NAME ${BOYLE_CXX_TEST_NAME} COMMAND ${BOYLE_CXX_TEST_NAME})
endfunction()
