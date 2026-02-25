local _scriptdir = os.scriptdir()

package("matplotplusplus")
    set_homepage("https://alandefreitas.github.io/matplotplusplus/")
    set_description("A C++ Graphics Library for Data Visualization")
    set_license("MIT")

    add_urls("https://github.com/alandefreitas/matplotplusplus/archive/refs/tags/v$(version).tar.gz")
    add_versions("1.2.2", "c7434b4fea0d0cc3508fd7104fafbb2fa7c824b1d2ccc51c52eaee26fc55a9a0")

    add_deps("cmake")
    add_links("matplot", "nodesoup")

    on_install(function (package)
        os.vrunv("patch", {"-p1", "--input=" .. path.join(_scriptdir, "disable_warning_as_error.patch")})
        local configs = {
            "-DCMAKE_CXX_STANDARD=23",
            "-DCMAKE_COMPILE_WARNING_AS_ERROR=OFF",
            "-DMATPLOTPP_BUILD_TESTS=OFF",
            "-DMATPLOTPP_BUILD_EXAMPLES=OFF",
            "-DMATPLOTPP_BUILD_WITH_JPEG=OFF",
            "-DMATPLOTPP_BUILD_WITH_TIFF=OFF",
            "-DMATPLOTPP_BUILD_WITH_PNG=OFF",
            "-DMATPLOTPP_BUILD_WITH_BLAS=OFF",
            "-DMATPLOTPP_BUILD_WITH_FFTW=OFF",
            "-DMATPLOTPP_BUILD_WITH_OPENCV=OFF",
        }
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:is_debug() and "Debug" or "Release"))
        table.insert(configs, "-DMATPLOTPP_BUILD_SHARED_LIBS=" .. (package:config("shared") and "ON" or "OFF"))
        import("package.tools.cmake").install(package, configs)
        local nodesoup_src = path.join(package:installdir("lib"), "Matplot++", "libnodesoup.a")
        if os.isfile(nodesoup_src) then
            os.cp(nodesoup_src, package:installdir("lib"))
        end
    end)

    on_test(function (package)
        assert(package:has_cxxincludes("matplot/matplot.h", {configs = {languages = "c++17"}}))
    end)
package_end()

add_requires("matplotplusplus 1.2.2")
