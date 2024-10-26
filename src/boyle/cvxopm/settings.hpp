/**
 * @file settings.hpp
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

namespace boyle::cvxopm {

template <std::floating_point Scalar = double, std::integral Index = int>
struct Settings {
    // Note: If this struct is updated, ensure update_settings and validate_settings are also
    // updated

    // Linear algebra settings
    Index device = 0;        // device identifier; currently used for CUDA devices
    Index linsys_solver = 1; // linear system solver to use

    // Control settings
    bool allocate_solution{true}; // boolean; allocate solution in OSQPSolver during osqp_setup
    Index verbose = 1;            // boolean; write out progress
    Index profiler_level = 0;     // integer; level of detail for profiler annotations
    bool warm_starting{true};     // boolean; warm start
    Index scaling = 10;           // data scaling iterations; if 0, then disabled
    bool polishing{false};        // boolean; polish ADMM solution

    // ADMM parameters
    Scalar rho = 0.1;      // ADMM penalty parameter
    bool rho_is_vec{true}; // boolean; is rho scalar or vector?
    Scalar sigma = 1E-6;   // ADMM penalty parameter
    Scalar alpha = 1.6;    // ADMM relaxation parameter

    // CG settings
    Scalar cg_max_iter = 20;       // maximum number of CG iterations per solve
    Scalar cg_tol_reduction = 10;  // number of consecutive zero CG iterations before the tolerance
                                   // gets halved
    Scalar cg_tol_fraction = 0.15; // CG tolerance (fraction of ADMM residuals)
    Index cg_precond = 1;          // Preconditioner to use in the CG method

    // adaptive rho logic
    bool adaptive_rho{true};         // boolean, is rho step size adaptive?
    Index adaptive_rho_interval = 0; // number of iterations between rho adaptations; if 0, then it
                                     // is timing-based
    Scalar adaptive_rho_fraction =
        0.4; // time interval for adapting rho (fraction of the setup time)
    Scalar adaptive_rho_tolerance = 5.0; // tolerance X for adapting rho; new rho must be X times
                                         // larger or smaller than the current one to change it

    // termination parameters
    Index max_iter = 4000;          // maximum number of iterations
    Scalar eps_abs = 1E-3;          // absolute solution tolerance
    Scalar eps_rel = 1E-3;          // relative solution tolerance
    Scalar eps_prim_inf = 1E-4;     // primal infeasibility tolerance
    Scalar eps_dual_inf = 1E-4;     // dual infeasibility tolerance
    bool scaled_termination{false}; // boolean; use scaled termination criteria
    Index check_termination = 25; // integer, check termination interval; if 0, checking is disabled
    Scalar time_limit = 1e10;     // maximum time to solve the problem (seconds)

    // polishing parameters
    Scalar delta = 1E-6;          // regularization parameter for polishing
    Index polish_refine_iter = 3; // number of iterative refinement steps in polishing
};

} // namespace boyle::cvxopm

namespace boost::serialization {

template <std::floating_point Scalar, std::integral Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::cvxopm::Settings<Scalar, Index>& obj,
    [[maybe_unused]] const unsigned int version
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

} // namespace boost::serialization
