toolchain("windows-msvc-x64")
    set_kind("standalone")
    set_homepage("https://visualstudio.microsoft.com/")
    set_description("Windows MSVC x86_64")

    on_check(function (toolchain)
        return is_host("windows") and is_arch("x86_64", "x64")
    end)

    on_load(function (toolchain)
        import("lib.detect.find_tool")

        toolchain:set("toolset", "cc", "cl.exe")
        toolchain:set("toolset", "cxx", "cl.exe")
        toolchain:set("toolset", "ld", "link.exe")
        toolchain:set("toolset", "sh", "link.exe")
        toolchain:set("toolset", "ar", "lib.exe")

        local ccache = find_tool("ccache")
        if ccache then
            toolchain:set("toolset", "ccache", ccache.program)
        end

        local base_cflags = {
            "/arch:AVX2", "/EHsc", "/utf-8", "/GS",
            "/guard:cf", "/wd4244", "/wd4267", "/wd4834",
            "/permissive-",
        }

        toolchain:add("cflags", table.unpack(base_cflags))
        toolchain:add("cxflags", table.unpack(base_cflags))
        toolchain:add("defines", "_USE_MATH_DEFINES")
    end)
toolchain_end()
