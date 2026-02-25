CPMAddPackage(
  NAME doctest
  VERSION 2.5.2
  URL https://github.com/doctest/doctest/releases/download/v2.5.2/doctest-v2.5.2.tar.gz
  URL_HASH SHA256=03635c39d59844f1550b9ffa56e3cbc79bddfed557e8465cd3692628aadba70a
  OPTIONS
    "DOCTEST_WITH_MAIN_IN_STATIC_LIB ON"
    "DOCTEST_NO_INSTALL ON"
)

target_compile_definitions(doctest INTERFACE DOCTEST_CONFIG_SUPER_FAST_ASSERTS)
target_compile_options(doctest INTERFACE
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-c2y-extensions>"
)
