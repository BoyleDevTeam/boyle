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
