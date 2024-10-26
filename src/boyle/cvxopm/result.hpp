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
#include <vector>

#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/aligned_allocator.hpp"

namespace boyle::cvxopm {

template <
    std::floating_point Scalar, typename Alloc = ::boyle::common::AlignedAllocator<Scalar, 32>>
struct [[nodiscard]] Result final {
    std::vector<Scalar, Alloc> prim_vars;
    std::vector<Scalar, Alloc> dual_vars;
    std::vector<Scalar, Alloc> prim_inf_cert;
    std::vector<Scalar, Alloc> dual_inf_cert;
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
