/**
 * @file result.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-11-03
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>

#include "boyle/common/utils/aligned_memory_resource.hpp"
#include "boyle/math/dense/vectorx.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar, std::integral Index = int>
struct [[nodiscard]] Result final {
    using value_type = Scalar;
    using index_type = Index;
    using param_type = ::boyle::math::pmr::VectorX<value_type>;
    using allocator_type = typename param_type::allocator_type;

    [[using gnu: always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return prim_vars.get_allocator();
    }

    param_type prim_vars{::boyle::common::pmr::getAlignedMemoryResource(32)};
    param_type dual_vars{::boyle::common::pmr::getAlignedMemoryResource(32)};
    param_type prim_inf_cert{::boyle::common::pmr::getAlignedMemoryResource(32)};
    param_type dual_inf_cert{::boyle::common::pmr::getAlignedMemoryResource(32)};
};

} // namespace boyle::cvxopm

namespace boost::serialization {

template <std::floating_point Scalar, typename Alloc>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::cvxopm::Result<Scalar, Alloc>& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.prim_vars;
    archive & obj.dual_vars;
    archive & obj.prim_inf_cert;
    archive & obj.dual_inf_cert;
    return;
}

} // namespace boost::serialization
