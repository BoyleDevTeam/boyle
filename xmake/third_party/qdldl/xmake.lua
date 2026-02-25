package("qdldl")
    set_homepage("https://github.com/osqp/qdldl")
    set_description("A free LDL factorisation routine")
    set_license("Apache-2.0")

    add_urls("https://github.com/osqp/qdldl/archive/refs/tags/v$(version).tar.gz")
    add_versions("0.1.9", "7d1285b2db15cf2730dc83b3d16ed28412f558591108cca4f28d4438bf72ceed")

    add_deps("cmake")

    on_install(function (package)
        local configs = {
            "-DCMAKE_POLICY_DEFAULT_CMP0069=NEW",
            "-DQDLDL_LONG=OFF",
            "-DQDLDL_FLOAT=OFF",
            "-DQDLDL_BUILD_DEMO_EXE=OFF",
            "-DQDLDL_UNITTESTS=OFF",
        }
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:is_debug() and "Debug" or "Release"))
        table.insert(configs, "-DQDLDL_BUILD_STATIC_LIB=" .. (package:config("shared") and "OFF" or "ON"))
        table.insert(configs, "-DQDLDL_BUILD_SHARED_LIB=" .. (package:config("shared") and "ON" or "OFF"))
        import("package.tools.cmake").install(package, configs)
    end)

    on_test(function (package)
        assert(package:has_cincludes("qdldl/qdldl.h"))
    end)
package_end()

add_requires("qdldl 0.1.9")
