/**
 * @file qp_problem.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-10-28
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <iterator>
#include <numeric>
#include <type_traits>
#include <unordered_map>

#include "boost/serialization/access.hpp"

#include "common/utils/macros.hpp"
#include "math/sparse_matrix/coo_matrix.hpp"
#include "math/sparse_matrix/lil_matrix.hpp"
#include "math/utils.hpp"

namespace tiny_pnc {
namespace math {

template <typename Scalar = double, typename Index = int>
class [[nodiscard]] QpProblem final {
    static_assert(std::is_arithmetic_v<Scalar>, "The loaded type must has arithmetic operators.");
    static_assert(std::is_integral_v<Index>, "The loaded type must be a floating-point type.");

    friend class boost::serialization::access;
    friend class OsqpSolver;

  public:
    [[using gnu: always_inline]] QpProblem(std::size_t num_vars, std::size_t num_cons) noexcept
        : objective_matrix_(num_vars, num_vars), objective_vector_(num_vars, 0.0),
          constrain_matrix_(num_cons, num_vars),
          lower_bounds_(num_cons, std::numeric_limits<Scalar>::lowest()),
          upper_bounds_(num_cons, std::numeric_limits<Scalar>::max()) {}

    ENABLE_IMPLICIT_CONSTRUCTORS(QpProblem);
    ~QpProblem() noexcept = default;

    [[using gnu: always_inline]]
    void resize(std::size_t num_vars, std::size_t num_cons) noexcept {
        objective_matrix_.resize(num_vars, num_vars);
        objective_vector_.resize(num_vars, 0.0);
        constrain_matrix_.resize(num_cons, num_vars);
        lower_bounds_.resize(num_cons, std::numeric_limits<Scalar>::lowest());
        upper_bounds_.resize(num_cons, std::numeric_limits<Scalar>::max());
        return;
    }

    [[using gnu: always_inline]]
    void clear() noexcept {
        objective_matrix_.clear();
        std::fill(objective_vector_.begin(), objective_vector_.end(), 0.0);
        constrain_matrix_.clear();
        std::fill(
            lower_bounds_.begin(), lower_bounds_.end(), std::numeric_limits<Scalar>::lowest()
        );
        std::fill(upper_bounds_.begin(), upper_bounds_.end(), std::numeric_limits<Scalar>::max());
        return;
    }

    [[using gnu: pure, always_inline]]
    std::size_t num_variables() const noexcept {
        return objective_matrix_.nrows();
    }

    [[using gnu: pure, always_inline]]
    std::size_t num_constraints() const noexcept {
        return constrain_matrix_.nrows();
    }

    template <typename ForwardIt>
    [[using gnu: pure]]
    Scalar cost(ForwardIt first, ForwardIt last) const noexcept {
        static_assert(
            std::is_same_v<typename std::iterator_traits<ForwardIt>::value_type, Scalar>,
            "The value_type of the iterators shall be the same as Scalar."
        );
        const std::size_t num_vars{static_cast<std::size_t>(last - first)};
        if (last - first != static_cast<long int>(objective_matrix_.nrows())) {
            spdlog::warn(
                "Mismatch size detected! The size of state_vec and num_variables has to be "
                "identical: state_vec.size() = {0:d} while num_variables() = {1:d}.",
                num_vars, objective_matrix_.nrows()
            );
            return 0.0;
        }
        Scalar cost{0.0};
        const std::vector<Scalar>& values{objective_matrix_.values()};
        const std::vector<Index>& row_indices{objective_matrix_.rowIndices()};
        const std::vector<Index>& col_indices{objective_matrix_.colIndices()};
        const std::size_t values_size{values.size()};
        for (std::size_t i{0}; i < values_size; ++i) {
            const Index row = row_indices[i];
            const Index col = col_indices[i];
            if (row != col) {
                cost += values[i] * first[row] * first[col];
            } else {
                cost += values[i] * first[row] * first[col] * 0.5;
            }
        }
        for (std::size_t i{0}; i < num_vars; ++i) {
            cost += objective_vector_[i] * first[i];
        }
        return cost;
    }

    template <typename ForwardIt>
    [[using gnu: pure]]
    bool validate(ForwardIt first, ForwardIt last) const noexcept {
        static_assert(
            std::is_same_v<typename std::iterator_traits<ForwardIt>::value_type, Scalar>,
            "The value_type of the iterators shall be the same as Scalar."
        );
        const std::size_t num_vars{static_cast<std::size_t>(last - first)};
        if (last - first != static_cast<long int>(objective_matrix_.nrows())) {
            spdlog::warn(
                "Mismatch size detected! The size of state_vec and num_variables has to be "
                "identical: state_vec.size() = {0:d} while num_variables() = {1:d}.",
                num_vars, objective_matrix_.nrows()
            );
            return false;
        }
        const std::size_t num_cons = constrain_matrix_.nrows();
        for (Index i{0}; i < static_cast<Index>(num_cons); ++i) {
            Scalar inner_prod{0.0};
            const std::vector<Scalar>& row_values_ = constrain_matrix_.rowValues(i);
            const std::vector<Index>& row_indices = constrain_matrix_.rowIndices(i);
            for (Index j{0}; j < static_cast<Index>(num_vars); ++j) {
                inner_prod += row_values_[j] * first[row_indices[j]];
            }
            if (inner_prod < lower_bounds_[i] || inner_prod > upper_bounds_[i]) {
                return false;
            }
        }
        return true;
    }

    [[using gnu: always_inline, hot]]
    void addQuadCostTerm(Index row, Index col, Scalar coeff) noexcept {
        if (row >= static_cast<Index>(objective_matrix_.nrows()) ||
            col >= static_cast<Index>(objective_matrix_.ncols())) {
            spdlog::warn(
                "Out-of-range error detected! The index ({0:d},{1:d}) is not in the allowed range: "
                "the allowed range is ({2:d},{3:d}).",
                row, col, objective_matrix_.nrows(), objective_matrix_.ncols()
            );
            return;
        }
        if (coeff == 0.0) {
            return;
        }
        if (row < col) {
            coeff += objective_matrix_.coeff(row, col);
            objective_matrix_.updateCoeff(row, col, coeff);
        } else if (row == col) {
            coeff *= 2.0;
            coeff += objective_matrix_.coeff(row, col);
            objective_matrix_.updateCoeff(row, col, coeff);
        } else {
            spdlog::warn(
                "The objective matrix has to be a upper triangular sparse matrix which requires "
                "row_index < col_index. row_index: {0:d}, col_index: {1:d}. This term is updated "
                "to ({1:d}, {0:d}) instead.",
                row, col
            );
            coeff += objective_matrix_.coeff(col, row);
            objective_matrix_.updateCoeff(col, row, coeff);
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    void updateQuadCostTerm(Index row, Index col, Scalar coeff) noexcept {
        if (row >= static_cast<Index>(objective_matrix_.nrows()) ||
            col >= static_cast<Index>(objective_matrix_.ncols())) {
            spdlog::warn(
                "Out-of-range error detected! The index ({0:d},{1:d}) is not in the allowed range: "
                "the allowed range is ({2:d},{3:d}).",
                row, col, objective_matrix_.nrows(), objective_matrix_.ncols()
            );
            return;
        }
        if (row < col) {
            objective_matrix_.updateCoeff(row, col, coeff);
        } else if (row == col) {
            objective_matrix_.updateCoeff(row, col, coeff * 2.0);
        } else {
            spdlog::warn(
                "The objective matrix has to be a upper triangular sparse matrix which requires "
                "row_index < col_index. row_index: {0:d}, col_index: {1:d}. This term is updated "
                "to ({1:d}, {0:d}) instead.",
                row, col
            );
            objective_matrix_.updateCoeff(col, row, coeff);
        }
        return;
    }

    [[using gnu: always_inline, hot]]
    void addLinCostTerm(Index row, Scalar coeff) noexcept {
        if (row >= static_cast<Index>(objective_vector_.size())) {
            spdlog::warn(
                "Out-of-range error detected! The input row has to be in the allowed range: row = "
                "{0:d} while max_row = {1:d}).",
                row, objective_vector_.size()
            );
            return;
        }
        objective_vector_[row] += coeff;
        return;
    }

    [[using gnu: always_inline, hot]]
    void updateLinCostTerm(Index row, Scalar coeff) noexcept {
        if (row >= static_cast<Index>(objective_vector_.size())) {
            spdlog::warn(
                "Out-of-range error detected! The input row has to be in the allowed range: row = "
                "{0:d} while max_row = {1:d}).",
                row, objective_vector_.size()
            );
            return;
        }
        objective_vector_[row] = coeff;
        return;
    }

    [[using gnu: always_inline, hot]]
    void addClampCostTerm(
        std::unordered_map<Index, Scalar> constrain_vec, Scalar offset, Scalar linear_coeff,
        Scalar quadratic_coeff = 0.0
    ) noexcept {
        const std::size_t num_vars = objective_matrix_.nrows();
        const std::size_t num_cons = constrain_matrix_.nrows();
        objective_matrix_.resize(num_vars + 1, num_vars + 1);
        objective_matrix_.updateCoeff(num_vars, num_vars, quadratic_coeff);
        objective_vector_.push_back(linear_coeff);
        constrain_matrix_.resize(num_cons + 2, num_vars + 1);
        constrain_matrix_.updateRow(num_cons, {{num_vars, 1.0}});
        lower_bounds_.push_back(0.0);
        upper_bounds_.push_back(std::numeric_limits<Scalar>::max());
        constrain_vec.emplace(num_vars, -1.0);
        constrain_matrix_.updateRow(num_cons + 1, constrain_vec);
        lower_bounds_.push_back(std::numeric_limits<Scalar>::lowest());
        upper_bounds_.push_back(offset);
        return;
    }

    [[using gnu: always_inline, hot]]
    void addConstrainTerm(
        const std::unordered_map<Index, Scalar>& constrain_vec, Scalar lower_bound,
        Scalar upper_bound
    ) noexcept {
        const std::size_t num_cons = constrain_matrix_.nrows();
        constrain_matrix_.resize(num_cons + 1, constrain_matrix_.ncols());
        constrain_matrix_.updateRow(num_cons, constrain_vec);
        lower_bounds_.push_back(lower_bound);
        upper_bounds_.push_back(upper_bound);
        return;
    }

    [[using gnu: always_inline, hot]]
    void updateConstrainTerm(
        Index row, const std::unordered_map<Index, Scalar>& constrain_vec, Scalar lower_bound,
        Scalar upper_bound
    ) noexcept {
        if (row >= static_cast<Index>(constrain_matrix_.nrows())) {
            spdlog::warn(
                "Out-of-range error detected! The input row has to be in the allowed range: row = "
                "{0:d} while max_row = {1:d}).",
                row, constrain_matrix_.nrows()
            );
            return;
        }
        constrain_matrix_.updateRow(row, constrain_vec);
        lower_bounds_[row] = lower_bound;
        upper_bounds_[row] = upper_bound;
        return;
    }

  private:
    template <class Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& objective_matrix_;
        ar& objective_vector_;
        ar& constrain_matrix_;
        ar& lower_bounds_;
        ar& upper_bounds_;
        return;
    }

    CooMatrix<Scalar, Index> objective_matrix_;
    std::vector<Scalar> objective_vector_;
    LilMatrix<Scalar, Index> constrain_matrix_;
    std::vector<Scalar> lower_bounds_;
    std::vector<Scalar> upper_bounds_;
};

} // namespace math
} // namespace tiny_pnc
