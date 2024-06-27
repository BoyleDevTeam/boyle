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
#include <format>
#include <stdexcept>

#include "boost/serialization/access.hpp"
#include "boost/serialization/array.hpp"
#include "boost/serialization/vector.hpp"
#include "osqp.h"

#include "boyle/common/utils/exec_on_exit.hpp"
#include "boyle/common/utils/macros.hpp"
#include "boyle/math/qp_solvers/qp_problem.hpp"
#include "boyle/math/sparse_matrix/csc_matrix.hpp"

namespace boyle::math {

using OsqpSettings = OSQPSettings;
using OsqpInfo = OSQPInfo;

struct [[nodiscard]] OsqpResult final {
    std::vector<OSQPFloat> prim_vars;
    std::vector<OSQPFloat> prim_inf_cert;
    std::vector<OSQPFloat> dual_vars;
    std::vector<OSQPFloat> dual_inf_cert;
    OsqpInfo info;
};

class [[nodiscard]] OsqpSolver final {
    friend class boost::serialization::access;

  public:
    [[using gnu: always_inline]] OsqpSolver() noexcept { osqp_set_default_settings(&settings); }
    [[using gnu: always_inline]] explicit OsqpSolver(const OsqpSettings& c_settings) noexcept
        : settings{c_settings} {}
    DISABLE_COPY_AND_MOVE(OsqpSolver);
    ~OsqpSolver() noexcept = default;

    template <std::floating_point Scalar, std::integral Index>
    [[using gnu: pure, flatten, leaf]] [[nodiscard]]
    auto solve(
        const QpProblem<Scalar, Index>& qp_problem, const std::vector<Scalar>& prim_vars_0 = {},
        const std::vector<Scalar>& dual_vars_0 = {}
    ) const -> OsqpResult {
#if BOYLE_CHECK_PARAMS == 1
        if (!prim_vars_0.empty() && prim_vars_0.size() != qp_problem.num_variables()) {
            throw std::invalid_argument(std::format(
                "Osqp runtime error detected! Size of prim_vars_0_ and num_vars must be identical: "
                "prim_vars_0_.size() = {0:d} while num_vars = {1:d}",
                prim_vars_0.size(), qp_problem.num_variables()
            ));
        }
        if (!dual_vars_0.empty() && dual_vars_0.size() != qp_problem.num_constraints()) {
            throw std::invalid_argument(std::format(
                "Osqp runtime error detected! Size of dual_vars_0_ and num_cons must be identical: "
                "dual_vars_0_.size() = {0:d} while num_cons = {1:d}",
                dual_vars_0.size(), qp_problem.num_constraints()
            ));
        }
#endif
        OSQPInt exit_flag{0};
        auto* solver = static_cast<OSQPSolver*>(malloc(sizeof(OSQPSolver)));
        auto* P = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));
        auto* A = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));

        ::boyle::common::ExecOnExit exec_on_exit{[&solver, &P, &A]() -> void {
            osqp_cleanup(solver);
            free(P);
            free(A);
            return;
        }};

        const auto num_vars = static_cast<OSQPInt>(qp_problem.num_variables());
        const auto num_cons = static_cast<OSQPInt>(qp_problem.num_constraints());
        CscMatrix<OSQPFloat, OSQPInt> objective_matrix{qp_problem.m_objective_matrix};
        std::vector<OSQPFloat> objective_vector{qp_problem.m_objective_vector};
        CscMatrix<OSQPFloat, OSQPInt> constrain_matrix{qp_problem.m_constrain_matrix};
        std::vector<OSQPFloat> lower_bounds{qp_problem.m_lower_bounds};
        std::vector<OSQPFloat> upper_bounds{qp_problem.m_upper_bounds};

        csc_set_data(
            P, num_vars, num_vars, static_cast<OSQPInt>(objective_matrix.nnzs()),
            const_cast<OSQPFloat*>(objective_matrix.values().data()),
            const_cast<OSQPInt*>(objective_matrix.innerIndices().data()),
            const_cast<OSQPInt*>(objective_matrix.outerIndices().data())
        );
        csc_set_data(
            A, num_cons, num_vars, static_cast<OSQPInt>(constrain_matrix.nnzs()),
            const_cast<OSQPFloat*>(constrain_matrix.values().data()),
            const_cast<OSQPInt*>(constrain_matrix.innerIndices().data()),
            const_cast<OSQPInt*>(constrain_matrix.outerIndices().data())
        );

        exit_flag = osqp_setup(
            &solver, P, objective_vector.data(), A, lower_bounds.data(), upper_bounds.data(),
            num_cons, num_vars, &settings
        );
        if (exit_flag) {
            throw std::runtime_error(std::format(
                "Osqp runtime error detected! Not able to setup a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            ));
        }

        // Set warm start
        if (!prim_vars_0.empty() && !dual_vars_0.empty()) {
            exit_flag = osqp_warm_start(solver, prim_vars_0.data(), dual_vars_0.data());
            if (exit_flag) {
                throw std::runtime_error(std::format(
                    "Osqp runtime error detected! Not able to set warm start for a OSQPSolver: "
                    "exit_flag = {0:d}.",
                    exit_flag
                ));
            }
        }

        exit_flag = osqp_solve(solver);
        if (exit_flag) {
            ;
            throw std::runtime_error(std::format(
                "Osqp runtime error detected! Not able to execute a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            ));
        }

        return OsqpResult{
            .prim_vars{solver->solution->x, solver->solution->x + num_vars},
            .prim_inf_cert{
                solver->solution->prim_inf_cert, solver->solution->prim_inf_cert + num_vars
            },
            .dual_vars{solver->solution->y, solver->solution->y + num_cons},
            .dual_inf_cert{
                solver->solution->dual_inf_cert, solver->solution->dual_inf_cert + num_cons
            },
            .info{*(solver->info)}
        };
    }

    OsqpSettings settings{};

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & settings;
        return;
    }
};

} // namespace boyle::math

namespace boost::serialization {

[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::math::OsqpSettings& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.device;
    archive & obj.linsys_solver;
    archive & obj.verbose;
    archive & obj.warm_starting;
    archive & obj.scaling;
    archive & obj.polishing;
    archive & obj.rho;
    archive & obj.rho_is_vec;
    archive & obj.sigma;
    archive & obj.alpha;
    archive & obj.cg_max_iter;
    archive & obj.cg_tol_reduction;
    archive & obj.cg_tol_fraction;
    archive & obj.cg_precond;
    archive & obj.adaptive_rho;
    archive & obj.adaptive_rho_interval;
    archive & obj.adaptive_rho_fraction;
    archive & obj.adaptive_rho_tolerance;
    archive & obj.max_iter;
    archive & obj.eps_abs;
    archive & obj.eps_rel;
    archive & obj.eps_prim_inf;
    archive & obj.eps_dual_inf;
    archive & obj.scaled_termination;
    archive & obj.check_termination;
    archive & obj.time_limit;
    archive & obj.delta;
    archive & obj.polish_refine_iter;
    return;
}

[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::math::OsqpInfo& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive& boost::serialization::make_array(obj.status, 32);
    archive & obj.status_val;
    archive & obj.status_polish;
    archive & obj.obj_val;
    archive & obj.prim_res;
    archive & obj.dual_res;
    archive & obj.iter;
    archive & obj.rho_updates;
    archive & obj.rho_estimate;
    archive & obj.setup_time;
    archive & obj.solve_time;
    archive & obj.update_time;
    archive & obj.polish_time;
    archive & obj.run_time;
    return;
}

[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::math::OsqpResult& obj, [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.prim_vars;
    archive & obj.prim_inf_cert;
    archive & obj.dual_vars;
    archive & obj.dual_inf_cert;
    archive & obj.info;
    return;
}

} // namespace boost::serialization
