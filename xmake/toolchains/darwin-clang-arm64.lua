toolchain("darwin-clang-arm64")
    set_kind("standalone")
    set_homepage("https://llvm.org/")
    set_description("Darwin Clang aarch64-apple-darwin")

    on_check(function (toolchain)
        return is_host("macosx") and is_arch("arm64")
    end)

    on_load(function (toolchain)
        import("lib.detect.find_tool")

        toolchain:set("toolset", "cc", "clang")
        toolchain:set("toolset", "cxx", "clang++")
        toolchain:set("toolset", "ld", "clang++")
        toolchain:set("toolset", "sh", "clang++")
        toolchain:set("toolset", "ar", "ar")

        local ccache = find_tool("ccache")
        if ccache then
            toolchain:set("toolset", "ccache", ccache.program)
        end

        local base_cflags = {
            "-march=native", "-pipe", "-fexceptions",
            "-fno-omit-frame-pointer", "-mno-omit-leaf-frame-pointer",
            "-Wall", "-Wextra", "-Wpedantic",
        }
        local base_cxxflags = {"-stdlib=libc++"}

        local llvm_prefix = try { function ()
            return os.iorun("brew --prefix llvm@21"):trim()
        end} or "/opt/homebrew/opt/llvm"

        local base_ldflags = {
            "-L" .. llvm_prefix .. "/lib/c++",
            "-L" .. llvm_prefix .. "/lib/unwind",
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
