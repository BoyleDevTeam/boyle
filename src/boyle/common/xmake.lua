includes("utils/xmake.lua")

target("common_fsm")
    set_kind("headeronly")
    add_headerfiles("fsm.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
target_end()
