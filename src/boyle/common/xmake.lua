includes("utils/xmake.lua")

target("common_fsm")
    set_kind("headeronly")
    add_includedirs("$(projectdir)/src", {public = true})
    add_headerfiles("fsm.hpp")
target_end()
