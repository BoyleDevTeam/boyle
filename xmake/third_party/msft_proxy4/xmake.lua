package("msft_proxy4")
    set_kind("library", {headeronly = true})
    set_homepage("https://github.com/microsoft/proxy")
    set_description("Proxy: Easy Polymorphism in C++")
    set_license("MIT")

    add_urls("https://github.com/ngcpp/proxy/releases/download/$(version)/proxy-$(version).tgz")
    add_versions("4.0.2", "6e24ceaf9b6dfb65c9e0bdab7dd5cb27dc88959d8d210d5f81a6842143735370")

    on_install(function (package)
        os.cp("proxy", package:installdir("include"))
    end)

    on_test(function (package)
        assert(package:has_cxxincludes("proxy/proxy.h", {configs = {languages = "c++20"}}))
    end)
package_end()

add_requires("msft_proxy4 4.0.2")
