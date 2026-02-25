toolchain("linux-clang-x64")
    set_kind("standalone")
    set_homepage("https://llvm.org/")
    set_description("Linux Clang x86_64-pc-linux-gnu")

    on_check(function (toolchain)
        return is_host("linux") and is_arch("x86_64")
    end)

    on_load(function (toolchain)
        import("lib.detect.find_tool")

        toolchain:set("toolset", "cc", "clang")
        toolchain:set("toolset", "cxx", "clang++")
        toolchain:set("toolset", "ld", "clang++")
        toolchain:set("toolset", "sh", "clang++")
        toolchain:set("toolset", "ar", "llvm-ar")
        toolchain:set("toolset", "ranlib", "llvm-ranlib")
        toolchain:set("toolset", "as", "clang")

        local ccache = find_tool("ccache")
        if ccache then
            toolchain:set("toolset", "ccache", ccache.program)
        end

        local lld = find_tool("ld.lld")
        if lld then
            toolchain:add("ldflags", "-fuse-ld=lld")
            toolchain:add("shflags", "-fuse-ld=lld")
        end

        local base_cflags = {
            "-march=native", "-pipe", "-fexceptions",
            "-fno-omit-frame-pointer", "-mno-omit-leaf-frame-pointer",
            "-Wall", "-Wextra", "-Wpedantic",
        }
        local base_cxxflags = {"-stdlib=libc++"}
        local base_ldflags = {
            "-stdlib=libc++", "-rtlib=compiler-rt",
            "--unwindlib=libunwind", "-lc++abi",
        }

        toolchain:add("cflags", table.unpack(base_cflags))
        toolchain:add("cxflags", table.unpack(base_cflags))
        toolchain:add("cxxflags", table.unpack(base_cxxflags))
        toolchain:add("ldflags", table.unpack(base_ldflags))
        toolchain:add("shflags", table.unpack(base_ldflags))
    end)
toolchain_end()
