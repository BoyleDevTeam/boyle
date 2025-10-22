/**
 * @file route_line_quintic_offset_model.hpp
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

#include <limits>
#include <span>
#include <vector>

#include "boyle/cvxopm/info.hpp"
#include "boyle/cvxopm/problems/qp_problem.hpp"
#include "boyle/cvxopm/settings.hpp"
#include "boyle/cvxopm/solvers/osqp_solver.hpp"
#include "boyle/kinetics/border2.hpp"
#include "boyle/kinetics/path2.hpp"
#include "boyle/math/curves/piecewise_quintic_curve.hpp"

namespace boyle::kinetics {

class [[nodiscard]] RouteLineQuinticOffsetModel final {
  public:
    RouteLineQuinticOffsetModel(const RouteLineQuinticOffsetModel& other) noexcept = delete;
    auto operator=(const RouteLineQuinticOffsetModel& other) noexcept
        -> RouteLineQuinticOffsetModel& = delete;
    RouteLineQuinticOffsetModel(RouteLineQuinticOffsetModel&& other) noexcept = delete;
    auto operator=(RouteLineQuinticOffsetModel&& other) noexcept
        -> RouteLineQuinticOffsetModel& = delete;
    ~RouteLineQuinticOffsetModel() noexcept = default;

    explicit RouteLineQuinticOffsetModel(
        std::vector<::boyle::math::Vec2d> sketch_points, std::vector<double> sample_ss
    ) noexcept(!BOYLE_CHECK_PARAMS);
    explicit RouteLineQuinticOffsetModel(
        std::vector<::boyle::math::Vec2d> sketch_points, std::vector<double> sample_ss,
        const ::boyle::cvxopm::Settings<double, int>& settings
    ) noexcept(!BOYLE_CHECK_PARAMS);
    auto num_samples() const noexcept -> std::size_t;
    auto qp_problem() const noexcept -> const ::boyle::cvxopm::QpProblem<double, int>&;
    auto settings() const noexcept -> const ::boyle::cvxopm::Settings<double, int>&;
    auto setHardBorders(std::span<const HardBorder2d> hard_borders) noexcept -> void;
    auto setSoftBorders(std::span<const SoftBorder2d> soft_borders) noexcept -> void;
    auto setDdxRange(double ddx_min, double ddx_max) noexcept -> void;
    auto setDdyRange(double ddy_min, double ddy_max) noexcept -> void;
    auto setD4xRange(double d4x_min, double d4x_max) noexcept -> void;
    auto setD4yRange(double d4y_min, double d4y_max) noexcept -> void;
    auto setInitialState(
        ::boyle::math::Vec2d r0,
        ::boyle::math::Vec2d t0 =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d n0 =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d j0 = {
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
        }
    ) noexcept -> void;
    auto setFinalState(
        ::boyle::math::Vec2d rf,
        ::boyle::math::Vec2d tf =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d nf =
            {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()},
        ::boyle::math::Vec2d jf = {
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()
        }
    ) noexcept -> void;
    auto setOffsetCost(double offset_weight) noexcept -> void;
    auto setCurvatureCost(double curvature_weight) noexcept -> void;
    auto setDCurvatureCost(double dcurature_weight) noexcept -> void;
    auto solve() const -> std::pair<Path2d, ::boyle::cvxopm::Info<double, int>>;
    auto clear() noexcept -> void;

  private:
    auto setIntegrationRelation() noexcept -> void;
    auto xIndex(int s_index) const noexcept -> int;
    auto ddxIndex(int s_index) const noexcept -> int;
    auto d4xIndex(int s_index) const noexcept -> int;
    auto yIndex(int s_index) const noexcept -> int;
    auto ddyIndex(int s_index) const noexcept -> int;
    auto d4yIndex(int s_index) const noexcept -> int;
    int m_num_samples{0};
    double m_length_scale{0.0};
    ::boyle::cvxopm::QpProblem<double, int> m_qp_problem{};
    ::boyle::cvxopm::Settings<double, int> m_settings{};
    ::boyle::math::PiecewiseQuinticCurve2d m_sketch_curve{};
    std::vector<::boyle::math::Vec2d> m_sample_points{};
    std::vector<double> m_sample_ss{};
    std::vector<double> m_hs{};
    std::vector<double> m_h2s{};
    std::vector<double> m_h3s{};
    std::vector<double> m_h4s{};
    std::vector<double> m_h6s{};
    std::vector<double> m_h8s{};
    std::vector<double> m_reciprocal_hs{};
};

} // namespace boyle::kinetics
