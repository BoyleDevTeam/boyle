function load(package)
    package:add("links", "openblas")
    package:add("defines", "HAVE_LAPACK_CONFIG_H", "LAPACK_COMPLEX_CPP")
end

function install(package)
    local configs = {
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        "-DCMAKE_COMPILE_WARNING_AS_ERROR=OFF",
        "-DBUILD_TESTING=OFF",
        "-DNOFORTRAN=ON",
        "-DC_LAPACK=ON",
        "-DDYNAMIC_ARCH=OFF",
        "-DUSE_THREAD=OFF",
        "-DUSE_OPENMP=OFF",
    }
    table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:is_debug() and "Debug" or "Release"))
    table.insert(configs, "-DBUILD_SHARED_LIBS=" .. (package:config("shared") and "ON" or "OFF"))
    import("package.tools.cmake").install(package, configs)
end
