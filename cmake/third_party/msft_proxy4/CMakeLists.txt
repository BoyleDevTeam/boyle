CPMAddPackage(
  NAME msft_proxy4
  VERSION 4.0.2
  URL https://github.com/ngcpp/proxy/releases/download/4.0.2/proxy-4.0.2.tgz
  URL_HASH SHA256=6e24ceaf9b6dfb65c9e0bdab7dd5cb27dc88959d8d210d5f81a6842143735370
  DOWNLOAD_ONLY True
)

add_library(msft_proxy4 INTERFACE)
target_sources(msft_proxy4
  INTERFACE
    FILE_SET HEADERS
      TYPE HEADERS
      BASE_DIRS ${msft_proxy4_SOURCE_DIR}
      FILES
        ${msft_proxy4_SOURCE_DIR}/proxy/proxy.h
        ${msft_proxy4_SOURCE_DIR}/proxy/proxy_fmt.h
        ${msft_proxy4_SOURCE_DIR}/proxy/proxy_macros.h
        ${msft_proxy4_SOURCE_DIR}/proxy/v4/proxy.h
        ${msft_proxy4_SOURCE_DIR}/proxy/v4/proxy_fmt.h
        ${msft_proxy4_SOURCE_DIR}/proxy/v4/proxy_macros.h
        ${msft_proxy4_SOURCE_DIR}/proxy/v4/proxy.ixx
)

set_target_properties(msft_proxy4
  PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${msft_proxy4_SOURCE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

if(BOYLE_ENABLE_INSTALL)
  install(
    TARGETS
      msft_proxy4
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()
