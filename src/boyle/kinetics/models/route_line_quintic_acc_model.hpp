/**
 * @file route_line_acc_model.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <vector>

#include "boyle/common/utils/macros.hpp"
#include "boyle/kinetics/fence1.hpp"
#include "boyle/kinetics/motion1.hpp"
#include "boyle/math/qp_solvers/osqp_solver.hpp"
#include "boyle/math/qp_solvers/qp_problem.hpp"

namespace boyle::kinetics {

class [[nodiscard]] RouteLineQuinticAccModel final {
  public:
    struct [[nodiscard]] Result final {
        Motion1d motion;
        ::boyle::math::OsqpInfo info;
    };

    explicit RouteLineQuinticAccModel(std::vector<double> sample_ts) noexcept(!BOYLE_CHECK_PARAMS);
    explicit RouteLineQuinticAccModel(
        std::vector<double> sample_ts, const ::boyle::math::OsqpSettings& settings
    ) noexcept(!BOYLE_CHECK_PARAMS);
    DISABLE_COPY_AND_MOVE(RouteLineQuinticAccModel);
    ~RouteLineQuinticAccModel() noexcept = default;
    auto num_samples() const noexcept -> std::size_t;
    auto qp_problem() const noexcept -> const ::boyle::math::QpProblem<double, int>&;
    auto settings() const noexcept -> const ::boyle::math::OsqpSettings&;
    auto setHardFences(const std::vector<HardFence1d>& hard_fences) noexcept -> void;
    auto setSoftFences(const std::vector<SoftFence1d>& soft_fences) noexcept -> void;
    auto setVelocityRange(double lower_bound, double upper_bound) noexcept -> void;
    auto setAccelRange(double lower_bound, double upper_bound) noexcept -> void;
    auto setInitialState(
        double s0 = std::numeric_limits<double>::quiet_NaN(),
        double v0 = std::numeric_limits<double>::quiet_NaN(),
        double a0 = std::numeric_limits<double>::quiet_NaN(),
        double j0 = std::numeric_limits<double>::quiet_NaN()
    ) noexcept -> void;
    auto setFinalState(
        double sf = std::numeric_limits<double>::quiet_NaN(),
        double vf = std::numeric_limits<double>::quiet_NaN(),
        double af = std::numeric_limits<double>::quiet_NaN(),
        double jf = std::numeric_limits<double>::quiet_NaN()
    ) noexcept -> void;
    auto setVelocityCost(double target_velocity, double velocity_weight) noexcept -> void;
    auto setAccelCost(double accel_weight) noexcept -> void;
    auto setJerkCost(double jerk_weight) noexcept -> void;
    auto setSnapCost(double snap_weight) noexcept -> void;
    auto solve() const -> Result;
    auto clear() noexcept -> void;

  private:
    auto setIntegrationRelation() noexcept -> void;
    auto sIndex(int t_index) const noexcept -> int;
    auto vIndex(int t_index) const noexcept -> int;
    auto aIndex(int t_index) const noexcept -> int;
    int m_num_samples{0};
    double m_time_scale{0.0};
    ::boyle::math::QpProblem<double, int> m_qp_problem{};
    ::boyle::math::OsqpSettings m_settings{};
    std::vector<double> m_sample_ts{};
    std::vector<double> m_hs{};
    std::vector<double> m_h2s{};
    std::vector<double> m_h3s{};
    std::vector<double> m_h4s{};
    std::vector<double> m_reciprocal_hs{};
    std::vector<double> m_reciprocal_h2s{};
    std::vector<double> m_reciprocal_h3s{};
    std::vector<double> m_reciprocal_h4s{};
};

} // namespace boyle::kinetics
