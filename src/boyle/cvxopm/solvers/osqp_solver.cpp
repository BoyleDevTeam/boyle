/**
 * @file osqp_solver.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-11-03
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/cvxopm/solvers/osqp_solver.hpp"

#include <format>
#include <memory>
#include <span>
#include <stdexcept>

#include "osqp.h"

#include "boyle/math/sparse/csc_matrix.hpp"

namespace boyle::cvxopm {

template <>
[[using gnu: pure]] [[nodiscard]]
auto OsqpSolver<OSQPFloat, OSQPInt>::solve(
    const QpProblem<OSQPFloat, OSQPInt>& qp_problem, std::span<const OSQPFloat> prim_vars_0,
    std::span<const OSQPFloat> dual_vars_0
) const -> std::pair<::boyle::cvxopm::Result<OSQPFloat>, Info<OSQPFloat, OSQPInt>> {
#if BOYLE_CHECK_PARAMS == 1
    if (!prim_vars_0.empty() && prim_vars_0.size() != qp_problem.num_variables()) [[unlikely]] {
        throw std::invalid_argument(
            std::format(
                "Osqp runtime error detected! Size of prim_vars_0_ and num_vars must be identical: "
                "prim_vars_0_.size() = {0:d} while num_vars = {1:d}",
                prim_vars_0.size(), qp_problem.num_variables()
            )
        );
    }
    if (!dual_vars_0.empty() && dual_vars_0.size() != qp_problem.num_constraints()) [[unlikely]] {
        throw std::invalid_argument(
            std::format(
                "Osqp runtime error detected! Size of dual_vars_0_ and num_cons must be identical: "
                "dual_vars_0_.size() = {0:d} while num_cons = {1:d}",
                dual_vars_0.size(), qp_problem.num_constraints()
            )
        );
    }
#endif

    OSQPSettings osqp_settings{
        .device = static_cast<OSQPInt>(settings.device),
        .linsys_solver = static_cast<enum osqp_linsys_solver_type>(settings.linsys_solver),
        .allocate_solution = static_cast<OSQPInt>(settings.allocate_solution),
        .verbose = static_cast<OSQPInt>(settings.verbose),
        .profiler_level = static_cast<OSQPInt>(settings.profiler_level),
        .warm_starting = static_cast<OSQPInt>(settings.warm_starting),
        .scaling = static_cast<OSQPInt>(settings.scaling),
        .polishing = static_cast<OSQPInt>(settings.polishing),
        .rho = static_cast<OSQPFloat>(settings.rho),
        .rho_is_vec = static_cast<OSQPInt>(settings.rho_is_vec),
        .sigma = static_cast<OSQPFloat>(settings.sigma),
        .alpha = static_cast<OSQPFloat>(settings.alpha),
        .cg_max_iter = static_cast<OSQPInt>(settings.cg_max_iter),
        .cg_tol_reduction = static_cast<OSQPInt>(settings.cg_tol_reduction),
        .cg_tol_fraction = static_cast<OSQPFloat>(settings.cg_tol_fraction),
        .cg_precond = static_cast<osqp_precond_type>(settings.cg_precond),
        .adaptive_rho = static_cast<OSQPInt>(settings.adaptive_rho),
        .adaptive_rho_interval = static_cast<OSQPInt>(settings.adaptive_rho_interval),
        .adaptive_rho_fraction = static_cast<OSQPFloat>(settings.adaptive_rho_fraction),
        .adaptive_rho_tolerance = static_cast<OSQPFloat>(settings.adaptive_rho_tolerance),
        .max_iter = static_cast<OSQPInt>(settings.max_iter),
        .eps_abs = static_cast<OSQPFloat>(settings.eps_abs),
        .eps_rel = static_cast<OSQPFloat>(settings.eps_rel),
        .eps_prim_inf = static_cast<OSQPFloat>(settings.eps_prim_inf),
        .eps_dual_inf = static_cast<OSQPFloat>(settings.eps_dual_inf),
        .scaled_termination = static_cast<OSQPInt>(settings.scaled_termination),
        .check_termination = static_cast<OSQPInt>(settings.check_termination),
        .check_dualgap = static_cast<OSQPInt>(settings.check_dualgap),
        .time_limit = static_cast<OSQPFloat>(settings.time_limit),
        .delta = static_cast<OSQPFloat>(settings.delta),
        .polish_refine_iter = static_cast<OSQPInt>(settings.polish_refine_iter)
    };

    const auto num_vars{static_cast<OSQPInt>(qp_problem.num_variables())};
    const auto num_cons{static_cast<OSQPInt>(qp_problem.num_constraints())};
    const ::boyle::math::pmr::CscMatrix<OSQPFloat, OSQPInt> objective_matrix{
        qp_problem.m_objective_matrix
    };
    const ::boyle::math::pmr::CscMatrix<OSQPFloat, OSQPInt> constrain_matrix{
        qp_problem.m_constrain_matrix
    };

    const OSQPCscMatrix P{
        .m = num_vars,
        .n = num_vars,
        .p = const_cast<OSQPInt*>(objective_matrix.outerIndices().data()),
        .i = const_cast<OSQPInt*>(objective_matrix.innerIndices().data()),
        .x = const_cast<OSQPFloat*>(objective_matrix.values().data()),
        .nzmax = static_cast<OSQPInt>(objective_matrix.nnzs()),
        .nz = -1,
        .owned = 0
    };

    const OSQPCscMatrix A{
        .m = num_cons,
        .n = num_vars,
        .p = const_cast<OSQPInt*>(constrain_matrix.outerIndices().data()),
        .i = const_cast<OSQPInt*>(constrain_matrix.innerIndices().data()),
        .x = const_cast<OSQPFloat*>(constrain_matrix.values().data()),
        .nzmax = static_cast<OSQPInt>(constrain_matrix.nnzs()),
        .nz = -1,
        .owned = 0
    };

    std::unique_ptr<OSQPSolver, decltype(&osqp_cleanup)> solver{nullptr, osqp_cleanup};

    OSQPInt exit_flag = osqp_setup(
        std::inout_ptr(solver), &P, qp_problem.m_objective_vector.data(), &A,
        qp_problem.m_lower_bounds.data(), qp_problem.m_upper_bounds.data(), num_cons, num_vars,
        &osqp_settings
    );

    if (exit_flag != 0) [[unlikely]] {
        throw std::runtime_error(
            std::format(
                "Osqp runtime error detected! Not able to setup a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            )
        );
    }

    // Set warm start
    if (!prim_vars_0.empty() && !dual_vars_0.empty()) {
        exit_flag = osqp_warm_start(solver.get(), prim_vars_0.data(), dual_vars_0.data());
        if (exit_flag != 0) [[unlikely]] {
            throw std::runtime_error(
                std::format(
                    "Osqp runtime error detected! Not able to set warm start for a OSQPSolver: "
                    "exit_flag = {0:d}.",
                    exit_flag
                )
            );
        }
    }

    exit_flag = osqp_solve(solver.get());
    if (exit_flag != 0) [[unlikely]] {
        throw std::runtime_error(
            std::format(
                "Osqp runtime error detected! Not able to execute a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            )
        );
    }

    ::boyle::cvxopm::Result<OSQPFloat> result{
        .prim_vars{solver->solution->x, solver->solution->x + num_vars, settings.memory_resource},
        .dual_vars{solver->solution->y, solver->solution->y + num_cons, settings.memory_resource},
        .prim_inf_cert{
            solver->solution->prim_inf_cert, solver->solution->prim_inf_cert + num_cons,
            settings.memory_resource
        },
        .dual_inf_cert{
            solver->solution->dual_inf_cert, solver->solution->dual_inf_cert + num_vars,
            settings.memory_resource
        },
    };

    ::boyle::cvxopm::Info<OSQPFloat, OSQPInt> info{
        .status{std::to_array(solver->info->status)},
        .status_val = solver->info->status_val,
        .status_polish = solver->info->status_polish,
        .obj_val = solver->info->obj_val,
        .dual_obj_val = solver->info->dual_obj_val,
        .prim_res = solver->info->prim_res,
        .dual_res = solver->info->dual_res,
        .duality_gap = solver->info->duality_gap,
        .iter = solver->info->iter,
        .rho_updates = solver->info->rho_updates,
        .rho_estimate = solver->info->rho_estimate,
        .setup_time = solver->info->setup_time,
        .solve_time = solver->info->solve_time,
        .update_time = solver->info->update_time,
        .polish_time = solver->info->polish_time,
        .run_time = solver->info->run_time,
        .primdual_int = solver->info->primdual_int,
        .rel_kkt_error = solver->info->rel_kkt_error,
    };

    return std::make_pair(std::move(result), info);
}

} // namespace boyle::cvxopm
