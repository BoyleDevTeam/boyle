after_build(function (target)
    local buildir = path.absolute("$(builddir)")
    local cc_json = path.join(buildir, "compile_commands.json")
    local link_path = path.absolute("compile_commands.json")
    if os.isfile(cc_json) and not os.isfile(link_path) then
        os.ln(cc_json, link_path)
    end
end)
