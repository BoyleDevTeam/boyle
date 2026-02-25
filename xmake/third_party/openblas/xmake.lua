local _scriptdir = os.scriptdir()

package("openblas")
    set_homepage("https://www.openblas.net/")
    set_description("An optimized BLAS library based on GotoBLAS2 1.13 BSD version")
    set_license("BSD-3-Clause")

    add_urls("https://github.com/OpenMathLib/OpenBLAS/releases/download/v$(version)/OpenBLAS-$(version).tar.gz")
    add_versions("0.3.32", "f8a1138e01fddca9e4c29f9684fd570ba39dedc9ca76055e1425d5d4b1a4a766")

    add_syslinks("pthread")
    add_includedirs("include", "include/openblas")

    on_load(function (package)
        if package:is_plat("linux") then
            import("linux", {rootdir = _scriptdir}).load(package)
        end
    end)

    on_install("linux", function (package)
        import("linux", {rootdir = _scriptdir}).install(package)
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
    add_requires("openblas 0.3.32")
end
