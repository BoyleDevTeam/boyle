toolchain("linux-gcc-x64")
    set_kind("standalone")
    set_homepage("https://gcc.gnu.org/")
    set_description("Linux GCC x86_64-pc-linux-gnu")

    on_check(function (toolchain)
        return is_host("linux") and is_arch("x86_64")
    end)

    on_load(function (toolchain)
        import("lib.detect.find_tool")

        toolchain:set("toolset", "cc", "gcc")
        toolchain:set("toolset", "cxx", "g++")
        toolchain:set("toolset", "ld", "g++")
        toolchain:set("toolset", "sh", "g++")
        toolchain:set("toolset", "ar", "ar")

        local ccache = find_tool("ccache")
        if ccache then
            toolchain:set("toolset", "ccache", ccache.program)
        end

        local bfd = find_tool("ld.bfd")
        if bfd then
            toolchain:add("ldflags", "-fuse-ld=bfd")
            toolchain:add("shflags", "-fuse-ld=bfd")
        end

        local base_cflags = {
            "-march=native", "-pipe", "-fno-plt", "-fexceptions",
            "-fstack-clash-protection", "-fcf-protection",
            "-fno-omit-frame-pointer", "-mno-omit-leaf-frame-pointer",
            "-Wall", "-Wextra", "-Wpedantic", "-Wno-maybe-uninitialized",
        }
        local base_cxxflags = {"-Wp,-D_GLIBCXX_ASSERTIONS"}

        toolchain:add("cflags", table.unpack(base_cflags))
        toolchain:add("cxflags", table.unpack(base_cflags))
        toolchain:add("cxxflags", table.unpack(base_cxxflags))
    end)
toolchain_end()
