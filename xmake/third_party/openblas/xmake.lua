local _scriptdir = os.scriptdir()

package("openblas")
    set_homepage("https://www.openblas.net/")
    set_description("An optimized BLAS library based on GotoBLAS2 1.13 BSD version")
    set_license("BSD-3-Clause")

    add_urls("https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(version)/OpenBLAS-$(version).tar.gz")
    add_versions("0.3.33", "6761af1d9f5d353ab4f0b7497be2643313b36c8f31caec0144bfef198e71e6ab")

    if is_plat("linux", "macosx") then
        add_syslinks("pthread")
    end
    add_includedirs("include", "include/openblas")

    on_load(function (package)
        if package:is_plat("linux") then
            import("linux", {rootdir = _scriptdir}).load(package)
        elseif package:is_plat("macosx") then
            import("darwin", {rootdir = _scriptdir}).load(package)
        elseif package:is_plat("windows") then
            import("windows", {rootdir = _scriptdir}).load(package)
        end
    end)

    on_install("linux", function (package)
        import("linux", {rootdir = _scriptdir}).install(package)
    end)

    on_install("macosx", function (package)
        import("darwin", {rootdir = _scriptdir}).install(package)
    end)

    on_install("windows", function (package)
        import("windows", {rootdir = _scriptdir}).install(package)
    end)

    on_test(function (package)
        assert(package:check_csnippets({test = [[
            void test() {
                double A[6] = {1.0, 2.0, 1.0, -3.0, 4.0, -1.0};
                double B[6] = {1.0, 2.0, 1.0, -3.0, 4.0, -1.0};
                double C[9] = {.5, .5, .5, .5, .5, .5, .5, .5, .5};
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 2, 1, A, 3, B, 3, 2, C, 3);
            }
        ]]}, {includes = "cblas.h"}))
    end)
package_end()

if is_config("boyle_use_blas_lapack", "OpenBLAS") then
    add_requires("openblas 0.3.33")
end

target("blas_lapack")
    set_kind("phony")
    if is_config("boyle_use_blas_lapack", "OpenBLAS") then
        add_packages("openblas", {public = true})
        add_defines("BOYLE_USE_BLAS_LAPACK", {public = true})
        if is_plat("windows") then
            add_defines("HAVE_LAPACK_CONFIG_H", "LAPACK_COMPLEX_CPP", {public = true})
            add_cxflags("/wd4190", {public = true})
        else
            add_cxflags("-Wno-c99-extensions", {public = true})
        end
    end
target_end()
