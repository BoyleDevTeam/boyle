package("doctest")
    set_kind("library", {headeronly = true})
    set_homepage("https://github.com/doctest/doctest")
    set_description("The fastest feature-rich C++11/14/17/20 single-header testing framework")
    set_license("MIT")

    add_urls("https://github.com/doctest/doctest/releases/download/v$(version)/doctest-v$(version).tar.gz")
    add_versions("2.5.2", "03635c39d59844f1550b9ffa56e3cbc79bddfed557e8465cd3692628aadba70a")

    add_deps("cmake")
    add_includedirs("include", "include/doctest")

    on_install(function (package)
        local configs = {
            "-DDOCTEST_WITH_TESTS=OFF",
            "-DDOCTEST_WITH_MAIN_IN_STATIC_LIB=OFF",
            "-DDOCTEST_NO_INSTALL=OFF",
        }
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:is_debug() and "Debug" or "Release"))
        import("package.tools.cmake").install(package, configs)
    end)

    on_test(function (package)
        assert(package:check_cxxsnippets({test = [[
            int factorial(int number) { return number <= 1 ? number : factorial(number - 1) * number; }
            TEST_CASE("testing the factorial function") {
                CHECK(factorial(1) == 1);
            }
        ]]}, {configs = {languages = "c++11"}, includes = "doctest.h", defines = "DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN"}))
    end)
package_end()

add_requires("doctest 2.5.2")
