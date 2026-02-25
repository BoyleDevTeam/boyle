/**
 * @file dense_problem.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2025-11-03
 *
 * @copyright Copyright (c) 2025 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>

#include "boyle/math/dense/vectorx.hpp"
#include "boyle/math/mdfunctions/mdfunction_proxy.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar, std::integral Index = int>
class DenseProblem final {
  public:
    using value_type = Scalar;
    using index_type = Index;
    using param_type = ::boyle::math::pmr::VectorX<value_type>;
    using size_type = typename param_type::size_type;

    DenseProblem() noexcept = default;
    DenseProblem(const DenseProblem& other) noexcept = default;
    auto operator=(const DenseProblem& other) noexcept -> DenseProblem& = default;
    DenseProblem(DenseProblem&& other) noexcept = default;
    auto operator=(DenseProblem&& other) noexcept -> DenseProblem& = default;
    ~DenseProblem() noexcept = default;

    explicit DenseProblem(::boyle::math::MdFunctionProxy<param_type> objective_function) noexcept
        : m_objective_function{std::move(objective_function)} {}

    [[using gnu: pure, always_inline]]
    auto cost(const param_type& x) const noexcept -> value_type {
        return m_objective_function->eval(x);
    }

    [[using gnu: pure, always_inline]]
    auto num_variables() const noexcept -> size_type {
        return m_objective_function->num_dimensions();
    }

    [[using gnu: pure, always_inline]]
    auto objective_function() const noexcept -> const ::boyle::math::MdFunctionProxy<param_type>& {
        return m_objective_function;
    }

  private:
    ::boyle::math::MdFunctionProxy<param_type> m_objective_function;
};

} // namespace boyle::cvxopm
