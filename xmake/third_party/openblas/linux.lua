function _detect_cc(package)
    local cc = package:build_getenv("cc")
    if cc then
        if cc:find("clang", 1, true) then
            return "clang", cc
        end
        if cc:find("gcc", 1, true) then
            return "gcc", cc
        end
    end
    return "gcc", cc or "cc"
end

function _detect_fc(cc_type)
    import("lib.detect.find_tool")
    local gfortran = find_tool("gfortran")
    local flang = find_tool("flang-new") or find_tool("flang")

    if cc_type == "gcc" then
        if gfortran then
            return "gfortran", gfortran.program
        end
    elseif cc_type == "clang" then
        if flang then
            return "llvmflang", flang.program
        end
        if gfortran then
            return "gfortran", gfortran.program
        end
    end

    if gfortran then
        return "gfortran", gfortran.program
    end
    if flang then
        return "llvmflang", flang.program
    end
    raise("openblas: no Fortran compiler found")
end

function _find_gfortran_runtime(package)
    local search_dirs = {"/usr/local/lib64", "/usr/lib"}
    for _, d in ipairs(os.dirs("/usr/lib/gcc/x86_64-linux-gnu/*") or {}) do
        table.insert(search_dirs, d)
    end

    local function find_lib(name)
        for _, dir in ipairs(search_dirs) do
            local p = path.join(dir, name)
            if os.isfile(p) then
                return p
            end
        end
        return nil
    end

    local gfortran_so = find_lib("libgfortran.so")
    if gfortran_so then
        os.cp(gfortran_so, package:installdir("lib"))
    end

    local quadmath_so = find_lib("libquadmath.so")
    if quadmath_so then
        os.cp(quadmath_so, package:installdir("lib"))
    end
end

function _find_flang_runtime(package)
    for _, llvm_dir in ipairs(os.dirs("/usr/lib/llvm-*") or {}) do
        for _, f in ipairs(os.files(path.join(llvm_dir, "lib/clang/*/lib/x86_64-pc-linux-gnu/libflang_rt.runtime.a")) or {}) do
            os.cp(f, package:installdir("lib"))
            break
        end
        local decimal = path.join(llvm_dir, "lib/libFortranDecimal.a")
        if os.isfile(decimal) then
            os.cp(decimal, package:installdir("lib"))
        end
    end
end

function load(package)
    local cc_type = _detect_cc(package)
    local fc_type = _detect_fc(cc_type)

    package:add("links", "openblas")
    if fc_type == "gfortran" then
        package:add("links", "gfortran", "quadmath")
    elseif fc_type == "llvmflang" then
        package:add("links", "flang_rt.runtime", "FortranDecimal")
    end
end

function install(package)
    local cc_type, cc = _detect_cc(package)
    local fc_type, fc = _detect_fc(cc_type)

    local cflags = ""
    local asmflags = ""
    local fflags = ""
    local ldflags = ""

    if cc_type == "gcc" and fc_type == "gfortran" then
        cflags = "-flto=auto -fno-fat-lto-objects"
        asmflags = cflags
        fflags = "-flto=auto -fno-fat-lto-objects"
        ldflags = "-fuse-ld=bfd"
    elseif cc_type == "clang" and fc_type == "llvmflang" then
        cflags = "-flto=thin"
        asmflags = cflags
        fflags = "-flto=thin"
        ldflags = "-fuse-ld=lld -rtlib=compiler-rt"
    elseif cc_type == "clang" and fc_type == "gfortran" then
        cflags = "-flto=thin"
        asmflags = cflags
        fflags = "-fno-lto"
        ldflags = "-fuse-ld=lld -rtlib=compiler-rt"
    end

    local make_configs = {
        "CC=" .. cc,
        "FC=" .. fc,
        "CFLAGS=" .. cflags,
        "ASMFLAGS=" .. asmflags,
        "FFLAGS=" .. fflags,
        "LDFLAGS=" .. ldflags,
        "DYNAMIC_ARCH=0",
        "NO_AFFINITY=0",
        "USE_THREAD=0",
        "USE_LOCKING=0",
        "USE_OPENMP=0",
        "NO_SHARED=1",
    }

    if cc_type == "clang" then
        local ar = package:build_getenv("ar")
        local ranlib = package:build_getenv("ranlib")
        if ar then
            table.insert(make_configs, "AR=" .. ar)
        end
        if ranlib then
            table.insert(make_configs, "RANLIB=" .. ranlib)
        end
    end

    import("package.tools.make")
    make.make(package, {"clean"})
    make.build(package, table.join({"libs", "netlib"}, make_configs))
    make.make(package, table.join({"PREFIX=" .. package:installdir(), "install"}, make_configs))

    if fc_type == "gfortran" then
        _find_gfortran_runtime(package)
    elseif fc_type == "llvmflang" then
        _find_flang_runtime(package)
    end
end
