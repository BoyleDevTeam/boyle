CPMAddPackage(
  NAME msft_proxy
  VERSION 4.0.1
  URL https://github.com/microsoft/proxy/releases/download/4.0.1/proxy-4.0.1.tgz
  URL_HASH SHA256=8b76d883782090e3e9871c96f89f6443242521ac867bef4be05672463dab1f8d
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
