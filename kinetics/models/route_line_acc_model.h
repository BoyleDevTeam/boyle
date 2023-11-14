/**
 * @file route_line_acc_model.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <numeric>
#include <vector>

#include "common/utils/macros.hpp"
#include "kinetics/fence1.hpp"
#include "kinetics/motion1.hpp"
#include "math/functions/piecewise_functions/piecewise_linear_function1.hpp"
#include "math/qp_solvers/osqp_solver.hpp"
#include "math/qp_solvers/qp_problem.hpp"

namespace tiny_pnc {
namespace kinetics {

class [[nodiscard]] RouteLineAccModel final {
  public:
    struct [[nodiscard]] Result {
        Motion1d motion;
        tiny_pnc::math::OsqpInfo info;
    };

    explicit RouteLineAccModel(std::vector<double> sample_ts);
    explicit RouteLineAccModel(
        std::vector<double> sample_ts, tiny_pnc::math::OsqpSettings settings
    );
    DISABLE_COPY_AND_MOVE(RouteLineAccModel);
    ~RouteLineAccModel() noexcept = default;
    const std::size_t& num_samples() const noexcept;
    void setSoftFences(const std::vector<SoftFence1d>& soft_fences) noexcept;
    void setHardFences(const std::vector<HardFence1d>& hard_fences) noexcept;
    void setVelocityRange(double lower_bound, double upper_bound) noexcept;
    void setAccelRange(double lower_bound, double upper_bound) noexcept;
    void setInitialState(
        double s0 = std::numeric_limits<double>::quiet_NaN(),
        double v0 = std::numeric_limits<double>::quiet_NaN(),
        double a0 = std::numeric_limits<double>::quiet_NaN()
    ) noexcept;
    void setFinalState(
        double sf = std::numeric_limits<double>::quiet_NaN(),
        double vf = std::numeric_limits<double>::quiet_NaN(),
        double af = std::numeric_limits<double>::quiet_NaN()
    ) noexcept;
    void setVelocityCost(double target_velocity, double velocity_weight) noexcept;
    void setAccelCost(double accel_weight) noexcept;
    void setJerkCost(double jerk_weight) noexcept;
    void setSnapCost(double snap_weight) noexcept;
    Result solve() const noexcept;

  private:
    static std::array<double, 6> prorationCoeffs(double ratio, double scale = 1.0);
    void setIntegrationRelation() noexcept;
    std::size_t sIndex(std::size_t t_index) const noexcept;
    std::size_t vIndex(std::size_t t_index) const noexcept;
    std::size_t aIndex(std::size_t t_index) const noexcept;
    std::size_t num_samples_;
    double time_scale_;
    tiny_pnc::math::QpProblem<double, int> qp_problem_;
    tiny_pnc::math::OsqpSettings settings_;
    std::vector<double> sample_ts_;
    std::vector<double> hs_;
    std::vector<double> h2s_;
    std::vector<double> h3s_;
    std::vector<double> h4s_;
    std::vector<double> reciprocal_hs_;
    std::vector<double> reciprocal_h2s_;
    std::vector<double> reciprocal_h3s_;
    std::vector<double> reciprocal_h4s_;
};

} // namespace kinetics
} // namespace tiny_pnc
