CPMAddPackage(
  NAME msft_proxy
  VERSION 4.0.0
  URL https://github.com/microsoft/proxy/releases/download/4.0.0/proxy-4.0.0.tgz
  URL_HASH SHA256=d49a35aea5e7f8a6ca2ea594909097d1c57c64e90b098fb129ca79cbe9220c16
  DOWNLOAD_ONLY True
)

add_library(msft_proxy INTERFACE)
target_sources(msft_proxy
  INTERFACE
    FILE_SET HEADERS
      TYPE HEADERS
      BASE_DIRS ${msft_proxy_SOURCE_DIR}
      FILES
        ${msft_proxy_SOURCE_DIR}/proxy/proxy.h
        ${msft_proxy_SOURCE_DIR}/proxy/proxy_fmt.h
        ${msft_proxy_SOURCE_DIR}/proxy/proxy_macros.h
        ${msft_proxy_SOURCE_DIR}/proxy/v4/proxy.h
        ${msft_proxy_SOURCE_DIR}/proxy/v4/proxy_fmt.h
        ${msft_proxy_SOURCE_DIR}/proxy/v4/proxy_macros.h
        ${msft_proxy_SOURCE_DIR}/proxy/v4/proxy.ixx
)

set_target_properties(msft_proxy
  PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${msft_proxy_SOURCE_DIR}>;$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

if(BOYLE_ENABLE_INSTALL)
  install(
    TARGETS
      msft_proxy
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()
