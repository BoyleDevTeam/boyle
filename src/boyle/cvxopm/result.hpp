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

namespace boyle::cvxopm {

template <std::floating_point Scalar = double, std::integral Index = int>
struct [[nodiscard]] Result final {
    std::vector<Scalar> prim_vars;
    std::vector<Scalar> prim_inf_cert;
    std::vector<Scalar> dual_vars;
    std::vector<Scalar> dual_inf_cert;
};

} // namespace boyle::cvxopm

namespace boost::serialization {

template <std::floating_point Scalar, std::integral Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::cvxopm::Result<Scalar, Index>& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.prim_vars;
    archive & obj.prim_inf_cert;
    archive & obj.dual_vars;
    archive & obj.dual_inf_cert;
    return;
}

} // namespace boost::serialization
