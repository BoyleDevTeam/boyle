add_requires("spdlog 1.17.0", {configs = {
    std_format = true,
    header_only = true,
}})

add_defines("SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE")
