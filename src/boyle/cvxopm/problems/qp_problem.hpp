/**
 * @file qp_problem.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-28
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <limits>
#include <span>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/unordered/unordered_flat_map.hpp"

#include "boyle/math/sparse/dok_matrix.hpp"
#include "boyle/math/sparse/lil_matrix.hpp"

namespace boyle::cvxopm {

template <std::floating_point Scalar = double, std::integral Index = int>
class [[nodiscard]] QpProblem final {
    friend class boost::serialization::access;
    template <std::floating_point, std::integral>
    friend class OsqpSolver;

  public:
    QpProblem() noexcept = default;
    QpProblem(const QpProblem& other) noexcept = default;
    auto operator=(const QpProblem& other) noexcept -> QpProblem& = default;
    QpProblem(QpProblem&& other) noexcept = default;
    auto operator=(QpProblem&& other) noexcept -> QpProblem& = default;
    ~QpProblem() noexcept = default;

    [[using gnu: always_inline]]
    QpProblem(std::size_t num_vars, std::size_t num_cons) noexcept
        : m_num_vars{static_cast<Index>(num_vars)}, m_num_cons{static_cast<Index>(num_cons)},
          m_objective_matrix{num_vars, num_vars}, m_objective_vector(num_vars, 0.0),
          m_constrain_matrix{num_cons, num_vars},
          m_lower_bounds(num_cons, std::numeric_limits<Scalar>::lowest()),
          m_upper_bounds(num_cons, std::numeric_limits<Scalar>::max()) {}

    [[using gnu: always_inline]]
    auto resize(std::size_t num_vars, std::size_t num_cons) noexcept -> void {
        m_num_vars = static_cast<Index>(num_vars);
        m_num_cons = static_cast<Index>(num_cons);
        m_objective_matrix.resize(num_vars, num_vars);
        m_objective_vector.resize(num_vars, 0.0);
        m_constrain_matrix.resize(num_cons, num_vars);
        m_lower_bounds.resize(num_cons, std::numeric_limits<Scalar>::lowest());
        m_upper_bounds.resize(num_cons, std::numeric_limits<Scalar>::max());
        return;
    }

    [[using gnu: always_inline]]
    auto clear() noexcept -> void {
        m_num_vars = 0;
        m_num_cons = 0;
        m_objective_matrix.clear();
        std::fill(m_objective_vector.begin(), m_objective_vector.end(), 0.0);
        m_constrain_matrix.clear();
        std::fill(
            m_lower_bounds.begin(), m_lower_bounds.end(), std::numeric_limits<Scalar>::lowest()
        );
        std::fill(m_upper_bounds.begin(), m_upper_bounds.end(), std::numeric_limits<Scalar>::max());
        return;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto num_variables() const noexcept -> std::size_t {
        return m_objective_matrix.nrows();
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    auto num_constraints() const noexcept -> std::size_t {
        return m_constrain_matrix.nrows();
    }

    [[using gnu: pure]]
    auto cost(std::span<const Scalar> x) const noexcept -> Scalar {
        const Index x_size = x.size();
        if (x_size != m_num_vars) [[unlikely]] {
            return 0.0;
        }
        Scalar cost{0.0};
        for (const auto& [index_pair, value] : m_objective_matrix.dictionary()) {
            const auto& [row, col] = index_pair;
            if (row != col) {
                cost += value * x[row] * x[col];
            } else {
                cost += value * x[row] * x[col] * 0.5;
            }
        }
        for (Index i{0}; i < m_num_vars; ++i) {
            cost += m_objective_vector[i] * x[i];
        }
        return cost;
    }

    [[using gnu: pure]]
    auto validate(std::span<const Scalar> x) const noexcept -> bool {
        const Index x_size = x.size();
        if (x_size != m_num_vars) [[unlikely]] {
            return false;
        }
        for (Index i{0}; i < m_num_cons; ++i) {
            Scalar inner_prod{0.0};
            for (const auto& [col, value] : m_constrain_matrix.row_dictionaries().at(i)) {
                inner_prod += value * x[col];
            }
            if (inner_prod < m_lower_bounds[i] || inner_prod > m_upper_bounds[i]) {
                return false;
            }
        }
        return true;
    }

    [[using gnu: always_inline, hot]]
    auto addQuadCostTerm(Index row, Index col, Scalar coeff) noexcept -> void {
        if (row >= m_num_vars || col >= m_num_vars) [[unlikely]] {
            return;
        }
        if (coeff == 0.0) [[unlikely]] {
            return;
        }
        if (row < col) {
            coeff += m_objective_matrix.coeff(row, col);
            m_objective_matrix.updateCoeff(row, col, coeff);
        } else if (row == col) {
            coeff *= 2.0;
            coeff += m_objective_matrix.coeff(row, col);
            m_objective_matrix.updateCoeff(row, col, coeff);
        } else {
            coeff += m_objective_matrix.coeff(col, row);
            m_objective_matrix.updateCoeff(col, row, coeff);
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    auto updateQuadCostTerm(Index row, Index col, Scalar coeff) noexcept -> void {
        if (row >= m_num_vars || col >= m_num_vars) [[unlikely]] {
            return;
        }
        if (row < col) {
            m_objective_matrix.updateCoeff(row, col, coeff);
        } else if (row == col) {
            m_objective_matrix.updateCoeff(row, col, coeff * 2.0);
        } else {
            m_objective_matrix.updateCoeff(col, row, coeff);
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    auto addLinCostTerm(Index row, Scalar coeff) noexcept -> void {
        if (row >= m_num_vars) [[unlikely]] {
            return;
        }
        m_objective_vector[row] += coeff;
        return;
    }

    [[using gnu: always_inline, hot]]
    auto updateLinCostTerm(Index row, Scalar coeff) noexcept -> void {
        if (row >= m_num_vars) [[unlikely]] {
            return;
        }
        m_objective_vector[row] = coeff;
        return;
    }

    [[using gnu: always_inline, hot]]
    auto addClampCostTerm(
        boost::unordered_flat_map<Index, Scalar> constrain_vec, Scalar offset, Scalar linear_coeff,
        Scalar quadratic_coeff = 0.0
    ) noexcept -> void {
        m_objective_matrix.resize(m_num_vars + 1, m_num_vars + 1);
        m_objective_matrix.updateCoeff(m_num_vars, m_num_vars, quadratic_coeff * 2.0);
        m_objective_vector.push_back(linear_coeff);
        m_constrain_matrix.resize(m_num_cons + 2, m_num_vars + 1);
        m_constrain_matrix.updateRow(m_num_cons, {{m_num_vars, 1.0}});
        m_lower_bounds.push_back(0.0);
        m_upper_bounds.push_back(std::numeric_limits<Scalar>::max());
        constrain_vec.emplace(m_num_vars, -1.0);
        m_constrain_matrix.updateRow(m_num_cons + 1, constrain_vec);
        m_lower_bounds.push_back(std::numeric_limits<Scalar>::lowest());
        m_upper_bounds.push_back(offset);
        m_num_vars += 1;
        m_num_cons += 2;
        return;
    }

    [[using gnu: always_inline, hot]]
    auto addConstrainTerm(
        const boost::unordered_flat_map<Index, Scalar>& constrain_vec, Scalar lower_bound,
        Scalar upper_bound
    ) noexcept -> void {
        m_constrain_matrix.resize(m_num_cons + 1, m_constrain_matrix.ncols());
        m_constrain_matrix.updateRow(m_num_cons, constrain_vec);
        m_lower_bounds.push_back(lower_bound);
        m_upper_bounds.push_back(upper_bound);
        m_num_cons += 1;
        return;
    }

    [[using gnu: always_inline, hot]]
    auto updateConstrainTerm(
        Index row, const boost::unordered_flat_map<Index, Scalar>& constrain_vec,
        Scalar lower_bound, Scalar upper_bound
    ) noexcept -> void {
        if (row >= m_num_cons) [[unlikely]] {
            return;
        }
        m_constrain_matrix.updateRow(row, constrain_vec);
        m_lower_bounds[row] = lower_bound;
        m_upper_bounds[row] = upper_bound;
        return;
    }

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_num_vars;
        archive & m_num_cons;
        archive & m_objective_matrix;
        archive & m_objective_vector;
        archive & m_constrain_matrix;
        archive & m_lower_bounds;
        archive & m_upper_bounds;
        return;
    }

    Index m_num_vars{0};
    Index m_num_cons{0};
    ::boyle::math::DokMatrix<Scalar, Index> m_objective_matrix{};
    std::vector<Scalar> m_objective_vector{};
    ::boyle::math::LilMatrix<Scalar, Index> m_constrain_matrix{};
    std::vector<Scalar> m_lower_bounds{};
    std::vector<Scalar> m_upper_bounds{};
};

} // namespace boyle::cvxopm
