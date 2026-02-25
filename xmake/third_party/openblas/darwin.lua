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
    return "clang", cc or "cc"
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

function _find_first(...)
    for _, p in ipairs({...}) do
        for _, f in ipairs(os.files(p) or {}) do
            return f
        end
    end
    return nil
end

function _find_gfortran_runtime(package)
    local gfortran = _find_first(
        "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.dylib",
        "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.a"
    )
    if gfortran then
        os.cp(gfortran, package:installdir("lib"))
    end

    local quadmath = _find_first(
        "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.dylib",
        "/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.a"
    )
    if quadmath then
        os.cp(quadmath, package:installdir("lib"))
    end
end

function _find_flang_runtime(package)
    for _, llvm_dir in ipairs(os.dirs("/opt/homebrew/opt/flang*") or {}) do
        for _, f in ipairs(os.files(path.join(llvm_dir, "lib/clang/*/lib/darwin/libflang_rt.runtime.dylib")) or {}) do
            os.cp(f, package:installdir("lib"))
            break
        end
        local decimal = path.join(llvm_dir, "lib/libFortranDecimal.dylib")
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

    if cc_type == "clang" then
        ldflags = "-rtlib=compiler-rt"
    end

    import("lib.detect.find_tool")
    local ranlib = try { function () return os.iorun("xcrun --find ranlib"):trim() end } or "ranlib"
    local ar = try { function () return os.iorun("xcrun --find ar"):trim() end } or "ar"

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
        "RANLIB=" .. ranlib,
        "AR=" .. ar,
    }

    import("package.tools.make")
    make.make(package, {"clean"})
    make.build(package, table.join("libs", make_configs))
    make.build(package, table.join("netlib", make_configs))
    make.make(package, table.join({"PREFIX=" .. package:installdir(), "install"}, make_configs))

    if fc_type == "gfortran" then
        _find_gfortran_runtime(package)
    elseif fc_type == "llvmflang" then
        _find_flang_runtime(package)
    end
end
