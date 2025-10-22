/**
 * @file info.hpp
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

#include <array>
#include <concepts>

#include "boost/serialization/array.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar, std::integral Index>
struct [[nodiscard]] Info final {
    // solver status
    std::array<char, 32> status; // Status string, e.g. 'solved'
    Index status_val;            // Status as Index.
    Index status_polish; // Polishing status: successful (1), unperformed (0), unsuccessful (-1)

    // solution quality
    Scalar obj_val;      // Primal objective v
    Scalar dual_obj_val; ///< Dual objective value
    Scalar prim_res;     // Norm of primal residual
    Scalar dual_res;     // Norm of dual residual
    Scalar duality_gap;  ///< Duality gap (Primal obj - Dual obj)

    // algorithm information
    Index iter;          // Number of iterations taken
    Index rho_updates;   // Number of rho updates performned
    Scalar rho_estimate; // Best rho estimate so far from residuals

    // timing information
    Scalar setup_time;  // Setup phase time (seconds)
    Scalar solve_time;  // Solve phase time (seconds)
    Scalar update_time; // Update phase time (seconds)
    Scalar polish_time; // Polish phase time (seconds)
    Scalar run_time;    // Total solve time (seconds)

    // Convergence information
    Scalar primdual_int;  ///< Integral of duality gap over time (Primal-dual integral), requires
                          ///< profiling
    Scalar rel_kkt_error; ///< Relative KKT error
};

} // namespace boyle::cvxopm

namespace boost::serialization {

template <std::floating_point Scalar, std::integral Index>
[[using gnu: always_inline]]
inline auto serialize(
    auto& archive, ::boyle::cvxopm::Info<Scalar, Index>& obj,
    [[maybe_unused]] const unsigned int version
) noexcept -> void {
    archive & obj.status;
    archive & obj.status_val;
    archive & obj.status_polish;
    archive & obj.obj_val;
    archive & obj.dual_obj_val;
    archive & obj.prim_res;
    archive & obj.dual_res;
    archive & obj.duality_gap;
    archive & obj.iter;
    archive & obj.rho_updates;
    archive & obj.rho_estimate;
    archive & obj.setup_time;
    archive & obj.solve_time;
    archive & obj.update_time;
    archive & obj.polish_time;
    archive & obj.run_time;
    archive & obj.primdual_int;
    archive & obj.rel_kkt_error;
    return;
}

} // namespace boost::serialization
