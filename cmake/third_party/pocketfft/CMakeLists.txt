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

target_include_directories(pocketfft INTERFACE
  "$<BUILD_INTERFACE:${pocketfft_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_compile_definitions(pocketfft INTERFACE POCKETFFT_NO_MULTITHREADING)

if(BOYLE_ENABLE_INSTALL)
  install(
    TARGETS
      pocketfft
    EXPORT ${PROJECT_NAME}Targets
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()
