/**
 * @file osqp_solver.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-29
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <span>

#include "boyle/cvxopm/info.hpp"
#include "boyle/cvxopm/problems/qp_problem.hpp"
#include "boyle/cvxopm/result.hpp"
#include "boyle/cvxopm/settings.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar, std::integral Index = int>
struct OsqpSolver final {
    OsqpSolver() noexcept = default;
    OsqpSolver(const OsqpSolver& other) noexcept = delete;
    auto operator=(const OsqpSolver& other) noexcept -> OsqpSolver& = delete;
    OsqpSolver(OsqpSolver&& other) noexcept = delete;
    auto operator=(OsqpSolver&& other) noexcept -> OsqpSolver& = delete;
    ~OsqpSolver() noexcept = default;

    [[using gnu: always_inline]]
    explicit OsqpSolver(const ::boyle::cvxopm::Settings<Scalar, Index>& c_settings) noexcept
        : settings{c_settings} {}

    [[using gnu: pure]] [[nodiscard]]
    auto solve(
        const QpProblem<Scalar, Index>& qp_problem, std::span<const Scalar> prim_vars_0 = {},
        std::span<const Scalar> dual_vars_0 = {}
    ) const -> std::pair<::boyle::cvxopm::Result<Scalar>, ::boyle::cvxopm::Info<Scalar, Index>>;

    ::boyle::cvxopm::Settings<Scalar, Index> settings;
};

} // namespace boyle::cvxopm

namespace boost::serialization {

template <std::floating_point Scalar, std::integral Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::cvxopm::OsqpSolver<Scalar, Index>& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.settings;
    return;
}

} // namespace boost::serialization
