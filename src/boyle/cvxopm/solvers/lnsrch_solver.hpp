/**
 * @file lnsrch_solver.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-24
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <iostream>
#include <utility>

#include "boyle/common/utils/aligned_memory_resource.hpp"
#include "boyle/cvxopm/info.hpp"
#include "boyle/cvxopm/problems/dense_problem.hpp"
#include "boyle/cvxopm/result.hpp"
#include "boyle/cvxopm/settings.hpp"
#include "boyle/cvxopm/solvers/detail/lnsrch.hpp"
#include "boyle/math/dense/vectorx.hpp"
#include "boyle/math/mdfunctions/mdfunction_proxy.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar, std::integral Index = int>
struct LnsrchSolver final {
    using value_type = Scalar;
    using index_type = Index;
    using param_type = ::boyle::math::pmr::VectorX<value_type>;
    using size_type = std::size_t;

    auto solve(
        const DenseProblem<value_type, index_type>& problem,
        param_type prim_vars_0 = param_type{::boyle::common::pmr::getAlignedMemoryResource(32)},
        [[maybe_unused]] param_type dual_vars_0 = param_type{
            ::boyle::common::pmr::getAlignedMemoryResource(32)
        }
    ) const -> std::pair<Result<value_type>, Info<value_type, index_type>> {
        const ::boyle::math::MdFunctionProxy<param_type>& objective_function{
            problem.objective_function()
        };
        const size_type num_vars{objective_function->num_dimensions()};
        const param_type zeros(num_vars, 0.0, settings.memory_resource.get());

        param_type x{zeros};
        if (!prim_vars_0.empty()) {
            std::copy_n(prim_vars_0.data(), num_vars, x.data());
        }
        param_type dx(
            num_vars, std::numeric_limits<value_type>::max(), settings.memory_resource.get()
        );
        index_type i{0};
        for (; i < settings.max_iter && !dx.identicalTo(zeros, settings.eps_abs); ++i, x += dx) {
            dx = detail::lnsrch<param_type>(
                objective_function, x, -objective_function->gradient(x), settings.max_step,
                settings.wolfe_rate, settings.eps_abs
            );
            // std::cout << dx[0] << " " << dx[1] << std::endl;
        }
        Result<value_type, index_type> result{.prim_vars{std::move(x)}};
        Info<value_type, index_type> info;
        info.iter = i;

        return std::make_pair(std::move(result), info);
    }

    Settings<value_type, index_type> settings;
};

} // namespace boyle::cvxopm
