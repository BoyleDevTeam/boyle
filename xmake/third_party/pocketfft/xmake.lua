local _scriptdir = os.scriptdir()

package("pocketfft")
    set_kind("library", {headeronly = true})
    set_homepage("https://github.com/mreineck/pocketfft")
    set_description("FFT implementation based on FFTPack, but with several improvements")
    set_license("BSD-3-Clause")

    add_urls("https://github.com/mreineck/pocketfft.git")
    add_versions("2023.02.14", "076cb3d2536b7c5d0629093ad886e10ac05f3623")

    add_patches("2023.02.14", path.join(_scriptdir, "add_cstdint.patch"), "d615b50c5e2961a17020fb7b94d9b9b9afe7ea5ddeee1fc2df5eccf5c439f1f1")

    add_defines("POCKETFFT_NO_MULTITHREADING")

    on_install(function (package)
        os.cp("pocketfft_hdronly.h", package:installdir("include"))
    end)

    on_test(function (package)
        assert(package:check_cxxsnippets({test = [[
            #include <pocketfft_hdronly.h>
            void test() {
                pocketfft::shape_t var;
            }
        ]]}, {configs = {languages = "c++11"}}))
    end)
package_end()

add_requires("pocketfft 2023.02.14")
