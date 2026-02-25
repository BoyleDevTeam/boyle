local _scriptdir = os.scriptdir()

package("doctest")
    set_kind("library", {headeronly = true})
    set_homepage("https://github.com/doctest/doctest")
    set_description("The fastest feature-rich C++11/14/17/20 single-header testing framework")
    set_license("MIT")

    add_urls("https://github.com/doctest/doctest/archive/refs/tags/v$(version).tar.gz")
    add_versions("2.4.12", "73381c7aa4dee704bd935609668cf41880ea7f19fa0504a200e13b74999c2d70")

    add_deps("cmake")
    add_includedirs("include", "include/doctest")

    on_install(function (package)
        os.vrunv("patch", {"-p1", "--input=" .. path.join(_scriptdir, "corrections_for_clang.patch")})
        local configs = {
            "-DCMAKE_POLICY_VERSION_MINIMUM=3.10",
            "-DCMAKE_COMPILE_WARNING_AS_ERROR=OFF",
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

add_requires("doctest 2.4.12")
