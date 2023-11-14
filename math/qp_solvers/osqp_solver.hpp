/**
 * @file osqp_solver.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-29
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

extern "C" {
#include "osqp.h"
}

#include "boost/serialization/access.hpp"
#include "boost/serialization/array.hpp"
#include "boost/serialization/vector.hpp"

#include "common/utils/exec_on_exit.hpp"
#include "common/utils/macros.hpp"
#include "math/qp_solvers/qp_problem.hpp"
#include "math/sparse_matrix/csc_matrix.hpp"

namespace tiny_pnc {
namespace math {

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

    [[using gnu: always_inline]] explicit OsqpSolver(const OsqpSettings& c_settings) noexcept {
        settings = c_settings;
    }

    DISABLE_COPY_AND_MOVE(OsqpSolver);

    ~OsqpSolver() noexcept = default;

    template <typename Scalar, typename Index>
    [[using gnu: pure, flatten, leaf]] [[nodiscard]]
    OsqpResult solve(
        QpProblem<Scalar, Index> qp_problem, const std::vector<Scalar>& prim_vars_0 = {},
        const std::vector<Scalar>& dual_vars_0 = {}
    ) const {
        if (!prim_vars_0.empty() && prim_vars_0.size() != qp_problem.num_variables()) {
            std::string error_msg = std::format(
                "Osqp runtime error detected! Size of prim_vars_0_ and num_vars must be identical: "
                "prim_vars_0_.size() = {0:d} while num_vars = {1:d}",
                prim_vars_0.size(), qp_problem.num_variables()
            );
            throw std::invalid_argument(std::move(error_msg));
        }
        if (!dual_vars_0.empty() && dual_vars_0.size() != qp_problem.num_constraints()) {
            std::string error_msg = std::format(
                "Osqp runtime error detected! Size of dual_vars_0_ and num_cons must be identical: "
                "dual_vars_0_.size() = {0:d} while num_cons = {1:d}",
                dual_vars_0.size(), qp_problem.num_constraints()
            );
            throw std::invalid_argument(std::move(error_msg));
        }

        OSQPSolver* solver = nullptr;
        OSQPCscMatrix* P = nullptr;
        OSQPCscMatrix* A = nullptr;
        OSQPInt exit_flag;
        const OSQPInt num_vars = static_cast<OSQPInt>(qp_problem.num_variables());
        const OSQPInt num_cons = static_cast<OSQPInt>(qp_problem.num_constraints());

        tiny_pnc::common::ExecOnExit exec_on_exit{[&solver, &P, &A]() -> void {
            osqp_cleanup(solver);
            free(P);
            free(A);
            return;
        }};

        solver = static_cast<OSQPSolver*>(malloc(sizeof(OSQPSolver)));
        P = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));
        A = static_cast<OSQPCscMatrix*>(malloc(sizeof(OSQPCscMatrix)));

        CscMatrix<OSQPFloat, OSQPInt> objective_matrix{qp_problem.objective_matrix_};
        CscMatrix<OSQPFloat, OSQPInt> constrain_matrix{qp_problem.constrain_matrix_};

        csc_set_data(
            P, num_vars, num_vars, objective_matrix.nnzs(),
            const_cast<OSQPFloat*>(objective_matrix.values().data()),
            const_cast<OSQPInt*>(objective_matrix.innerIndices().data()),
            const_cast<OSQPInt*>(objective_matrix.outerIndices().data())
        );
        csc_set_data(
            A, num_cons, num_vars, constrain_matrix.nnzs(),
            const_cast<OSQPFloat*>(constrain_matrix.values().data()),
            const_cast<OSQPInt*>(constrain_matrix.innerIndices().data()),
            const_cast<OSQPInt*>(constrain_matrix.outerIndices().data())
        );

        exit_flag = osqp_setup(
            &solver, P, qp_problem.objective_vector_.data(), A, qp_problem.lower_bounds_.data(),
            qp_problem.upper_bounds_.data(), num_cons, num_vars, &settings
        );
        if (exit_flag) {
            std::string error_msg = std::format(
                "Osqp runtime error detected! Not able to setup a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            );
            throw std::runtime_error(std::move(error_msg));
        }

        // Set warm start
        if (!prim_vars_0.empty() && !dual_vars_0.empty()) {
            exit_flag = osqp_warm_start(solver, prim_vars_0.data(), dual_vars_0.data());
            if (exit_flag) {
                std::string error_msg = std::format(
                    "Osqp runtime error detected! Not able to set warm start for a OSQPSolver: "
                    "exit_flag = {0:d}.",
                    exit_flag
                );
                throw std::runtime_error(std::move(error_msg));
            }
        }

        exit_flag = osqp_solve(solver);
        if (exit_flag) {
            std::string error_msg = std::format(
                "Osqp runtime error detected! Not able to execute a OSQPSolver: exit_flag = {0:d}.",
                exit_flag
            );
            throw std::runtime_error(std::move(error_msg));
        }

        return OsqpResult{
            .prim_vars{solver->solution->x, solver->solution->x + num_vars},
            .prim_inf_cert{
                solver->solution->prim_inf_cert, solver->solution->prim_inf_cert + num_vars},
            .dual_vars{solver->solution->y, solver->solution->y + num_cons},
            .dual_inf_cert{
                solver->solution->dual_inf_cert, solver->solution->dual_inf_cert + num_cons},
            .info{*(solver->info)}};
    }

    OsqpSettings settings;

  private:
    template <typename Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& settings;
        return;
    }
};

} // namespace math
} // namespace tiny_pnc

namespace boost {
namespace serialization {

template <typename Archive>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::math::OsqpSettings& obj, const unsigned int version
) noexcept {
    ar& obj.device;
    ar& obj.linsys_solver;
    ar& obj.verbose;
    ar& obj.warm_starting;
    ar& obj.scaling;
    ar& obj.polishing;
    ar& obj.rho;
    ar& obj.rho_is_vec;
    ar& obj.sigma;
    ar& obj.alpha;
    ar& obj.cg_max_iter;
    ar& obj.cg_tol_reduction;
    ar& obj.cg_tol_fraction;
    ar& obj.cg_precond;
    ar& obj.adaptive_rho;
    ar& obj.adaptive_rho_interval;
    ar& obj.adaptive_rho_fraction;
    ar& obj.adaptive_rho_tolerance;
    ar& obj.max_iter;
    ar& obj.eps_abs;
    ar& obj.eps_rel;
    ar& obj.eps_prim_inf;
    ar& obj.eps_dual_inf;
    ar& obj.scaled_termination;
    ar& obj.check_termination;
    ar& obj.time_limit;
    ar& obj.delta;
    ar& obj.polish_refine_iter;
    return;
}

template <typename Archive>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::math::OsqpInfo& obj, const unsigned int version
) noexcept {
    ar& boost::serialization::make_array(obj.status, 32);
    ar& obj.status_val;
    ar& obj.status_polish;
    ar& obj.obj_val;
    ar& obj.prim_res;
    ar& obj.dual_res;
    ar& obj.iter;
    ar& obj.rho_updates;
    ar& obj.rho_estimate;
    ar& obj.setup_time;
    ar& obj.solve_time;
    ar& obj.update_time;
    ar& obj.polish_time;
    ar& obj.run_time;
    return;
}

template <typename Archive>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, tiny_pnc::math::OsqpResult& obj, const unsigned int version
) noexcept {
    ar& obj.prim_vars;
    ar& obj.prim_inf_cert;
    ar& obj.dual_vars;
    ar& obj.dual_inf_cert;
    ar& obj.info;
    return;
}

} // namespace serialization
} // namespace boost
