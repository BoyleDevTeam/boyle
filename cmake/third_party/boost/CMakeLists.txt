CPMAddPackage(
  NAME Boost
  VERSION 1.91.0-1
  URL https://github.com/boostorg/boost/releases/download/boost-1.91.0-1/boost-1.91.0-1-cmake.tar.xz
  URL_HASH SHA256=cc5dc5006ecbdf0051f90979be31b4eee5987d9ae14ae9fb9c03cfa43fa3cdad
  OPTIONS
    "CMAKE_COMPILE_WARNING_AS_ERROR OFF"
    "CMAKE_UNITY_BUILD OFF"
    "BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS}"
    "BOOST_RUNTIME_LINK $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
    "BOOST_ENABLE_CMAKE ON"
    "BOOST_SKIP_INSTALL_RULES OFF"
    "BOOST_INCLUDE_LIBRARIES headers\\\;serialization"
)
