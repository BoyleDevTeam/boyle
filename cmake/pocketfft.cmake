CPMAddPackage(
  NAME pocketfft
  GITHUB_REPOSITORY mreineck/pocketfft
  GIT_TAG cpp
  DOWNLOAD_ONLY True
)

add_library(pocketfft INTERFACE)
target_sources(pocketfft
  INTERFACE
    FILE_SET HEADERS
      TYPE HEADERS
      BASE_DIRS ${pocketfft_SOURCE_DIR}
      FILES ${pocketfft_SOURCE_DIR}/pocketfft_hdronly.h
)

set_target_properties(pocketfft
  PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${pocketfft_SOURCE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    INTERFACE_COMPILE_DEFINITIONS "POCKETFFT_NO_MULTITHREADING"
)

if(BOYLE_ENABLE_INSTALL)
  install(
    TARGETS
      pocketfft
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()
