/**
 * @file route_line_cubic_offset_model.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-01-10
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <vector>

#include "boyle/common/utils/macros.hpp"
#include "boyle/cvxopm/problems/qp_problem.hpp"
#include "boyle/cvxopm/solvers/osqp_solver.hpp"
#include "boyle/kinetics/border2.hpp"
#include "boyle/kinetics/path2.hpp"
#include "boyle/math/curves/piecewise_curves/piecewise_quintic_curve2.hpp"

namespace boyle::kinetics {

class [[nodiscard]] RouteLineCubicOffsetModel final {
  public:
    explicit RouteLineCubicOffsetModel(
        const std::vector<::boyle::math::Vec2d>& sketch_points, std::vector<double> sample_ss
    ) noexcept(!BOYLE_CHECK_PARAMS);
    explicit RouteLineCubicOffsetModel(
        const std::vector<::boyle::math::Vec2d>& sketch_points, std::vector<double> sample_ss,
        const ::boyle::cvxopm::OsqpSolver::Settings& settings
    ) noexcept(!BOYLE_CHECK_PARAMS);
    DISABLE_COPY_AND_MOVE(RouteLineCubicOffsetModel);
    ~RouteLineCubicOffsetModel() noexcept = default;
    auto num_samples() const noexcept -> std::size_t;
    auto qp_problem() const noexcept -> const ::boyle::cvxopm::QpProblem<double, int>&;
    auto settings() const noexcept -> const ::boyle::cvxopm::OsqpSolver::Settings&;
    auto setHardBorders(const std::vector<HardBorder2d>& hard_borders) noexcept -> void;
    auto setSoftBorders(const std::vector<SoftBorder2d>& soft_borders) noexcept -> void;
    auto setDdxRange(double ddx_min, double ddx_max) noexcept -> void;
    auto setDdyRange(double ddy_min, double ddy_max) noexcept -> void;
    auto setInitialState(
        ::boyle::math::Vec2d r0,
        ::boyle::math::Vec2d t0 =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d n0 =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}
    ) noexcept -> void;
    auto setFinalState(
        ::boyle::math::Vec2d rf,
        ::boyle::math::Vec2d tf =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d nf =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}
    ) noexcept -> void;
    auto setOffsetCost(double offset_weight) noexcept -> void;
    auto setCurvatureCost(double curvature_weight) noexcept -> void;
    auto setDCurvatureCost(double dcurvature_weight) noexcept -> void;
    auto solve() const -> std::pair<Path2d, ::boyle::cvxopm::OsqpSolver::Info>;
    auto clear() noexcept -> void;

  private:
    auto setIntegrationRelation() noexcept -> void;
    auto xIndex(int s_index) const noexcept -> int;
    auto ddxIndex(int s_index) const noexcept -> int;
    auto yIndex(int s_index) const noexcept -> int;
    auto ddyIndex(int s_index) const noexcept -> int;
    int m_num_samples{0};
    double m_length_scale{0.0};
    ::boyle::cvxopm::QpProblem<double, int> m_qp_problem{};
    ::boyle::cvxopm::OsqpSolver::Settings m_settings{};
    ::boyle::math::PiecewiseQuinticCurve2d m_sketch_curve{};
    std::vector<::boyle::math::Vec2d> m_sample_points{};
    std::vector<double> m_sample_ss{};
    std::vector<double> m_hs{};
    std::vector<double> m_h2s{};
    std::vector<double> m_h4s{};
    std::vector<double> m_reciprocal_hs{};
};

} // namespace boyle::kinetics
