target("math_index_pair")
    set_kind("headeronly")
    add_headerfiles("index_pair.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
        add_packages("boost", {public = true})
    end
target_end()

target("math_sparse_traits")
    set_kind("headeronly")
    add_headerfiles("sparse_traits.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_deps("common_aligned_allocator", "math_concepts", "math_index_pair")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
        add_packages("boost", {public = true})
    end
target_end()

target("math_sparse_matrix_proxy")
    set_kind("headeronly")
    add_headerfiles("sparse_matrix_proxy.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("msft_proxy4", {public = true})
    add_deps("math_concepts")
target_end()

target("math_dok_matrix")
    set_kind("headeronly")
    add_headerfiles("dok_matrix.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("boost", {public = true})
    add_deps("math_concepts", "math_sparse_traits")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
    end
target_end()

target("math_lil_matrix")
    set_kind("headeronly")
    add_headerfiles("lil_matrix.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("boost", {public = true})
    add_deps("common_aligned_allocator", "math_concepts", "math_sparse_traits")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
    end
target_end()

target("math_coo_matrix")
    set_kind("headeronly")
    add_headerfiles("coo_matrix.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("boost", {public = true})
    add_deps("math_concepts", "math_sparse_traits")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
    end
target_end()

target("math_csc_matrix")
    set_kind("headeronly")
    add_headerfiles("csc_matrix.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("boost", {public = true})
    add_deps("math_concepts", "math_index_pair", "math_sparse_traits")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
    end
target_end()

target("math_csr_matrix")
    set_kind("headeronly")
    add_headerfiles("csr_matrix.hpp")
    add_includedirs("$(projectdir)/src", {public = true})
    add_packages("boost", {public = true})
    add_deps("math_concepts", "math_index_pair", "math_sparse_traits")
    if is_config("boyle_use_boost_unordered", true) then
        add_defines("BOYLE_USE_BOOST_UNORDERED=1", {public = true})
    end
target_end()
