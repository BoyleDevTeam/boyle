/**
 * @file route_line_acc_model.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/kinetics/models/route_line_quintic_acc_model.hpp"

#include <format>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "boyle/common/utils/logging.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_linear_function1.hpp"

namespace {

[[using gnu: const, always_inline, hot]] [[nodiscard]]
constexpr auto prorationCoeffs(double ratio, double scale) noexcept -> std::array<double, 6> {
    const double ratio2 = ratio * ratio;
    const double ratio3 = ratio2 * ratio;
    const double ratio4 = ratio3 * ratio;
    const double ratio5 = ratio4 * ratio;
    const double scale2 = scale * scale;
    return std::array<double, 6>{
        1 - ratio3 * 10.0 + ratio4 * 15.0 - ratio5 * 6.0,
        ratio3 * 10.0 - ratio4 * 15.0 + ratio5 * 6.0,
        (ratio - ratio3 * 6.0 + ratio4 * 8.0 - ratio5 * 3.0) * scale,
        -(ratio3 * 4.0 - ratio4 * 7.0 + ratio5 * 3.0) * scale,
        (ratio2 - ratio3 * 3.0 + ratio4 * 3.0 - ratio5) * scale2 * 0.5,
        (ratio3 - ratio4 * 2.0 + ratio5) * scale2 * 0.5
    };
}

} // namespace

namespace boyle::kinetics {

RouteLineQuinticAccModel::RouteLineQuinticAccModel(std::vector<double> sample_ts
) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
    if (sample_ts.size() < 2) {
        throw std::invalid_argument(std::format(
            "Invalid argument error detected: size of sample_ts must larger than 2: "
            "sample_ts.size() = {0:d}.",
            sample_ts.size()
        ));
    }
#endif
    m_num_samples = static_cast<int>(sample_ts.size());
    m_qp_problem.resize(m_num_samples * 3, m_num_samples * 5);
    osqp_set_default_settings(&m_settings);
    m_settings.scaling = 0;
    m_sample_ts = std::move(sample_ts);
    m_time_scale = m_sample_ts.back() - m_sample_ts.front();
    m_hs.reserve(m_num_samples - 1);
    m_h2s.reserve(m_num_samples - 1);
    m_h3s.reserve(m_num_samples - 1);
    m_h4s.reserve(m_num_samples - 1);
    m_reciprocal_hs.reserve(m_num_samples - 1);
    m_reciprocal_h2s.reserve(m_num_samples - 1);
    m_reciprocal_h3s.reserve(m_num_samples - 1);
    m_reciprocal_h4s.reserve(m_num_samples - 1);
    for (std::vector<double>::const_iterator it{m_sample_ts.cbegin() + 1}; it != m_sample_ts.cend();
         ++it) {
        const double diff = *it - *(it - 1);
        const double reciprocal_diff = 1.0 / diff;
        m_hs.push_back(diff);
        m_h2s.push_back(diff * diff);
        m_h3s.push_back(m_h2s.back() * diff);
        m_h4s.push_back(m_h3s.back() * diff);
        m_reciprocal_hs.push_back(reciprocal_diff);
        m_reciprocal_h2s.push_back(reciprocal_diff * reciprocal_diff);
        m_reciprocal_h3s.push_back(m_reciprocal_h2s.back() * reciprocal_diff);
        m_reciprocal_h4s.push_back(m_reciprocal_h3s.back() * reciprocal_diff);
    }
    setIntegrationRelation();
}

RouteLineQuinticAccModel::RouteLineQuinticAccModel(
    std::vector<double> sample_ts, const ::boyle::cvxopm::OsqpSolver::Settings& settings
) noexcept(!BOYLE_CHECK_PARAMS)
    : RouteLineQuinticAccModel{std::move(sample_ts)} {
    m_settings = settings;
}

auto RouteLineQuinticAccModel::num_samples() const noexcept -> std::size_t {
    return m_sample_ts.size();
}

auto RouteLineQuinticAccModel::qp_problem() const noexcept
    -> const ::boyle::cvxopm::QpProblem<double, int>& {
    return m_qp_problem;
}

auto RouteLineQuinticAccModel::settings() const noexcept
    -> const ::boyle::cvxopm::OsqpSolver::Settings& {
    return m_settings;
}

auto RouteLineQuinticAccModel::setIntegrationRelation() noexcept -> void {
    for (int i{1}; i < m_num_samples - 1; ++i) {
        std::unordered_map<int, double> constrain_vec{
            {sIndex(i - 1), m_reciprocal_h3s[i - 1] * 20.0},
            {sIndex(i), -(m_reciprocal_h3s[i - 1] + m_reciprocal_h3s[i]) * 20.0},
            {sIndex(i + 1), m_reciprocal_h3s[i] * 20.0},
            {vIndex(i - 1), m_reciprocal_h2s[i - 1] * 8.0},
            {vIndex(i), (m_reciprocal_h2s[i - 1] - m_reciprocal_h2s[i]) * 12.0},
            {vIndex(i + 1), -m_reciprocal_h2s[i] * 8.0},
            {aIndex(i - 1), m_reciprocal_hs[i - 1]},
            {aIndex(i), -(m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 3.0},
            {aIndex(i + 1), m_reciprocal_hs[i]}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 3 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
        constrain_vec = {
            {sIndex(i - 1), m_reciprocal_h4s[i - 1] * 30.0},
            {sIndex(i), -(m_reciprocal_h4s[i - 1] - m_reciprocal_h4s[i]) * 30.0},
            {sIndex(i + 1), -m_reciprocal_h4s[i] * 30.0},
            {vIndex(i - 1), m_reciprocal_h3s[i - 1] * 14.0},
            {vIndex(i), (m_reciprocal_h3s[i - 1] + m_reciprocal_h3s[i]) * 16.0},
            {vIndex(i + 1), m_reciprocal_h3s[i] * 14.0},
            {aIndex(i - 1), m_reciprocal_h2s[i - 1] * 2.0},
            {aIndex(i), -(m_reciprocal_h2s[i - 1] - m_reciprocal_h2s[i]) * 3.0},
            {aIndex(i + 1), -m_reciprocal_h2s[i] * 2.0}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 4 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
    }
    return;
}

auto RouteLineQuinticAccModel::setHardFences(const std::vector<HardFence1d>& hard_fences
) noexcept -> void {
    std::vector<double> lower_bound(m_num_samples, std::numeric_limits<double>::lowest());
    std::vector<double> upper_bound(m_num_samples, std::numeric_limits<double>::max());
    for (const HardFence1d& hard_fence : hard_fences) {
        const int istart = ::boyle::math::nearestUpperElement(
                               m_sample_ts.cbegin(), m_sample_ts.cend(), hard_fence.bound_ts.front()
                           ) -
                           m_sample_ts.cbegin();
        const int iend = ::boyle::math::nearestUpperElement(
                             m_sample_ts.cbegin(), m_sample_ts.cend(), hard_fence.bound_ts.back()
                         ) -
                         m_sample_ts.cbegin();
        if (istart == m_num_samples) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The front() of hard_fence.bound_ts should be "
                "less than m_sample_ts.back(): hard_fence.bound_ts.front() = {0:f} while "
                "m_sample_ts.back() = {1:f}.",
                hard_fence.bound_ts.front(), m_sample_ts.back()
            );
            break;
        }
        if (iend == 0) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The back() of hard_fence.bound_ts should be "
                "larger than m_sample_ts.front(): hard_fence.bound_ts.back() = {0:f} while "
                "m_sample_ts.front() = {1:f}.",
                hard_fence.bound_ts.back(), m_sample_ts.front()
            );
            break;
        }

        const ::boyle::math::PiecewiseLinearFunction1d bound_func{
            hard_fence.bound_ts, hard_fence.bound_ss
        };

        if (hard_fence.actio == ::boyle::kinetics::Actio::PUSHING) {
            if (istart != 0) {
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                const double ratio{(hard_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, hard_fence.bound_ss.front(), std::numeric_limits<double>::max()
                );
            }
            for (int i{istart}; i < iend; ++i) {
                lower_bound[i] = std::max(bound_func(m_sample_ts[i]), lower_bound[i]);
            }
            if (iend != m_num_samples) {
                const double h{m_sample_ts[iend] - m_sample_ts[iend - 1]};
                const double ratio{(hard_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, hard_fence.bound_ss.back(), std::numeric_limits<double>::max()
                );
            }
        } else {
            if (istart != 0) {
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                const double ratio{(hard_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(),
                    hard_fence.bound_ss.front()
                );
            }
            for (int i{istart}; i < iend; ++i) {
                upper_bound[i] = std::min(bound_func(m_sample_ts[i]), upper_bound[i]);
            }
            if (iend != m_num_samples) {
                const double h{m_sample_ts[iend] - m_sample_ts[iend - 1]};
                const double ratio{(hard_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(), hard_fence.bound_ss.back()
                );
            }
        }
    }

    for (int i{1}; i < m_num_samples - 1; ++i) {
        m_qp_problem.updateConstrainTerm(
            sIndex(i), {{sIndex(i), 1.0}}, lower_bound[i], upper_bound[i]
        );
    }
    return;
}

auto RouteLineQuinticAccModel::setSoftFences(const std::vector<SoftFence1d>& soft_fences
) noexcept -> void {
    for (const SoftFence1d& soft_fence : soft_fences) {
        const int istart = ::boyle::math::nearestUpperElement(
                               m_sample_ts.cbegin(), m_sample_ts.cend(), soft_fence.bound_ts.front()
                           ) -
                           m_sample_ts.cbegin();
        const int iend = ::boyle::math::nearestUpperElement(
                             m_sample_ts.cbegin(), m_sample_ts.cend(), soft_fence.bound_ts.back()
                         ) -
                         m_sample_ts.cbegin();
        if (istart == m_num_samples) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The front() of hard_fence.bound_ts should be "
                "less than m_sample_ts.back(): hard_fence.bound_ts.front() = {0:f} while "
                "m_sample_ts.back() = {1:f}.",
                soft_fence.bound_ts.front(), m_sample_ts.back()
            );
            break;
        }
        if (iend == 0) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! The back() of hard_fence.bound_ts should be "
                "larger than m_sample_ts.front(): hard_fence.bound_ts.back() = {0:f} while "
                "m_sample_ts.front() = {1:f}.",
                soft_fence.bound_ts.back(), m_sample_ts.front()
            );
            break;
        }

        const ::boyle::math::PiecewiseLinearFunction1d bound_func{
            soft_fence.bound_ts, soft_fence.bound_ss
        };

        double factor{0.0};

        if (soft_fence.actio == ::boyle::kinetics::Actio::PUSHING) {
            if (istart == iend) {
                factor =
                    (soft_fence.bound_ts.back() - soft_fence.bound_ts.front()) / m_time_scale * 0.5;
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                double ratio{(soft_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), -proration_coeffs[0]},
                    {sIndex(istart), -proration_coeffs[1]},
                    {vIndex(istart - 1), -proration_coeffs[2]},
                    {vIndex(istart), -proration_coeffs[3]},
                    {aIndex(istart - 1), -proration_coeffs[4]},
                    {aIndex(istart), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -soft_fence.bound_ss.front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                ratio = (soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h;
                proration_coeffs = prorationCoeffs(ratio, h);
                constrain_vec = {
                    {sIndex(iend - 1), -proration_coeffs[0]}, {sIndex(iend), -proration_coeffs[1]},
                    {vIndex(iend - 1), -proration_coeffs[2]}, {vIndex(iend), -proration_coeffs[3]},
                    {aIndex(iend - 1), -proration_coeffs[4]}, {aIndex(iend), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -soft_fence.bound_ss.back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                factor = (m_sample_ts[istart] - soft_fence.bound_ts.front()) / m_time_scale;
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                const double ratio{(soft_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), -proration_coeffs[0]},
                    {sIndex(istart), -proration_coeffs[1]},
                    {vIndex(istart - 1), -proration_coeffs[2]},
                    {vIndex(istart), -proration_coeffs[3]},
                    {aIndex(istart - 1), -proration_coeffs[4]},
                    {aIndex(istart), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -soft_fence.bound_ss.front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (int i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = m_hs[i] / m_time_scale * 0.5;
                } else if (i == iend - 1) {
                    factor = m_hs[i - 1] / m_time_scale * 0.5;
                } else {
                    factor = (m_hs[i - 1] + m_hs[i]) / m_time_scale * 0.5;
                }
                m_qp_problem.addClampCostTerm(
                    {{sIndex(i), -1.0}}, -bound_func(m_sample_ts[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != m_num_samples) {
                factor = (soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / m_time_scale;
                const double h{m_sample_ts[iend] - m_sample_ts[iend - 1]};
                const double ratio{(soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), -proration_coeffs[0]}, {sIndex(iend), -proration_coeffs[1]},
                    {vIndex(iend - 1), -proration_coeffs[2]}, {vIndex(iend), -proration_coeffs[3]},
                    {aIndex(iend - 1), -proration_coeffs[4]}, {aIndex(iend), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -soft_fence.bound_ss.back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        } else {
            if (istart == iend) {
                factor =
                    (soft_fence.bound_ts.back() - soft_fence.bound_ts.front()) / m_time_scale * 0.5;
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                double ratio{(soft_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), soft_fence.bound_ss.front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                ratio = (soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h;
                proration_coeffs = prorationCoeffs(ratio, h);
                constrain_vec = {
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), soft_fence.bound_ss.back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                factor = (m_sample_ts[istart] - soft_fence.bound_ts.front()) / m_time_scale;
                const double h{m_sample_ts[istart] - m_sample_ts[istart - 1]};
                const double ratio{(soft_fence.bound_ts.front() - m_sample_ts[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), soft_fence.bound_ss.front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (int i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = m_hs[i] / m_time_scale * 0.5;
                } else if (i == iend - 1) {
                    factor = m_hs[i - 1] / m_time_scale * 0.5;
                } else {
                    factor = (m_hs[i - 1] + m_hs[i]) / m_time_scale * 0.5;
                }
                m_qp_problem.addClampCostTerm(
                    {{sIndex(i), 1.0}}, bound_func(m_sample_ts[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != m_num_samples) {
                factor = (soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / m_time_scale;
                const double h{m_sample_ts[iend] - m_sample_ts[iend - 1]};
                const double ratio{(soft_fence.bound_ts.back() - m_sample_ts[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), soft_fence.bound_ss.back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        }
    }
    return;
}

auto RouteLineQuinticAccModel::setVelocityRange(double lower_bound, double upper_bound) noexcept
    -> void {
    for (int i{1}; i < m_num_samples - 1; ++i) {
        m_qp_problem.updateConstrainTerm(vIndex(i), {{vIndex(i), 1.0}}, lower_bound, upper_bound);
    }
    return;
}

auto RouteLineQuinticAccModel::setAccelRange(double lower_bound, double upper_bound) noexcept
    -> void {
    for (int i{1}; i < m_num_samples - 1; ++i) {
        m_qp_problem.updateConstrainTerm(aIndex(i), {{aIndex(i), 1.0}}, lower_bound, upper_bound);
    }
    return;
}

auto RouteLineQuinticAccModel::setInitialState(double s0, double v0, double a0, double j0) noexcept
    -> void {
    if (!std::isnan(s0)) {
        m_qp_problem.updateConstrainTerm(
            sIndex(0), {{sIndex(0), 1.0}}, s0 - ::boyle::math::kEpsilon,
            s0 + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(v0)) {
        m_qp_problem.updateConstrainTerm(
            vIndex(0), {{vIndex(0), 1.0}}, v0 - ::boyle::math::kEpsilon,
            v0 + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(a0)) {
        m_qp_problem.updateConstrainTerm(
            aIndex(0), {{aIndex(0), 1.0}}, a0 - ::boyle::math::kEpsilon,
            a0 + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(j0)) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 3,
            {{sIndex(0), -m_reciprocal_h3s[0] * 60.0},
             {sIndex(1), m_reciprocal_h3s[0] * 60.0},
             {vIndex(0), -m_reciprocal_h2s[0] * 36.0},
             {vIndex(1), -m_reciprocal_h2s[0] * 24.0},
             {aIndex(0), -m_reciprocal_hs[0] * 9.0},
             {aIndex(1), m_reciprocal_hs[0] * 3.0}},
            j0 - ::boyle::math::kEpsilon, j0 + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        m_num_samples * 4,
        {{sIndex(0), m_reciprocal_h4s[0] * 360.0},
         {sIndex(1), -m_reciprocal_h4s[0] * 360.0},
         {vIndex(0), m_reciprocal_h3s[0] * 192.0},
         {vIndex(1), m_reciprocal_h3s[0] * 168.0},
         {aIndex(0), m_reciprocal_h2s[0] * 36.0},
         {aIndex(1), -m_reciprocal_h2s[0] * 24.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );*/
    return;
}

auto RouteLineQuinticAccModel::setFinalState(double sf, double vf, double af, double jf) noexcept
    -> void {
    if (!std::isnan(sf)) {
        m_qp_problem.updateConstrainTerm(
            sIndex(m_num_samples - 1), {{sIndex(m_num_samples - 1), 1.0}},
            sf - ::boyle::math::kEpsilon, sf + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(vf)) {
        m_qp_problem.updateConstrainTerm(
            vIndex(m_num_samples - 1), {{vIndex(m_num_samples - 1), 1.0}},
            vf - ::boyle::math::kEpsilon, vf + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(af)) {
        m_qp_problem.updateConstrainTerm(
            aIndex(m_num_samples - 1), {{aIndex(m_num_samples - 1), 1.0}},
            af - ::boyle::math::kEpsilon, af + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(jf)) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 4 - 1,
            {{sIndex(m_num_samples - 2), -m_reciprocal_h3s[m_num_samples - 2] * 60.0},
             {sIndex(m_num_samples - 1), m_reciprocal_h3s[m_num_samples - 2] * 60.0},
             {vIndex(m_num_samples - 2), -m_reciprocal_h2s[m_num_samples - 2] * 24.0},
             {vIndex(m_num_samples - 1), -m_reciprocal_h2s[m_num_samples - 2] * 36.0},
             {aIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2] * 3.0},
             {aIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2] * 9.0}},
            jf - ::boyle::math::kEpsilon, jf + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        m_num_samples * 5 - 1,
        {{sIndex(m_num_samples - 2), -m_reciprocal_h4s[m_num_samples - 2] * 360.0},
         {sIndex(m_num_samples - 1), m_reciprocal_h4s[m_num_samples - 2] * 360.0},
         {vIndex(m_num_samples - 2), -m_reciprocal_h3s[m_num_samples - 2] * 168.0},
         {vIndex(m_num_samples - 1), -m_reciprocal_h3s[m_num_samples - 2] * 192.0},
         {aIndex(m_num_samples - 2), -m_reciprocal_h2s[m_num_samples - 2] * 24.0},
         {aIndex(m_num_samples - 1), m_reciprocal_h2s[m_num_samples - 2] * 36.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );*/
    return;
}

auto RouteLineQuinticAccModel::setVelocityCost(
    double target_velocity, double velocity_weight
) noexcept -> void {
    constexpr std::array<double, 9> kFactors{10.0 / 7.0, 20.0 / 7.0,  3.0 / 7.0,
                                             1.0 / 42.0, 8.0 / 35.0,  1.0 / 35.0,
                                             1.0 / 30.0, 1.0 / 105.0, 1.0 / 630.0};
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor = velocity_weight * m_reciprocal_hs[i] / m_time_scale;
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * kFactors[1]);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i), factor * m_hs[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * m_hs[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i), factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * m_hs[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * m_hs[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i), factor * m_h2s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i + 1), -factor * m_h2s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i), factor * m_h3s[i] * kFactors[6]);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i + 1), factor * m_h3s[i] * kFactors[7]);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * m_h2s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i), -factor * m_h3s[i] * kFactors[7]);
        m_qp_problem.addQuadCostTerm(
            vIndex(i + 1), aIndex(i + 1), -factor * m_h3s[i] * kFactors[6]
        );
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i), factor * m_h4s[i] * kFactors[8]);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i + 1), factor * m_h4s[i] * kFactors[8]);
        m_qp_problem.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * m_h4s[i] * kFactors[8]);
        m_qp_problem.addLinCostTerm(sIndex(i), factor * target_velocity * m_hs[i] * 2.0);
        m_qp_problem.addLinCostTerm(sIndex(i + 1), -factor * target_velocity * m_hs[i] * 2.0);
    }
    return;
}

auto RouteLineQuinticAccModel::setAccelCost(double accel_weight) noexcept -> void {
    constexpr std::array<double, 9> kFactors{120.0 / 7.0,  240.0 / 7.0,  6.0 / 7.0,
                                             192.0 / 35.0, 216.0 / 35.0, 22.0 / 35.0,
                                             8.0 / 35.0,   3.0 / 35.0,   1.0 / 35.0};
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor = accel_weight * m_reciprocal_h3s[i] / m_time_scale;
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * kFactors[1]);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i), factor * m_hs[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * m_hs[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i), factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * m_hs[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * m_hs[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i), factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * m_h2s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i), factor * m_h3s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * m_h3s[i] * kFactors[6]);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * m_h2s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * m_h3s[i] * kFactors[6]);
        m_qp_problem.addQuadCostTerm(
            vIndex(i + 1), aIndex(i + 1), -factor * m_h3s[i] * kFactors[5]
        );
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i), factor * m_h4s[i] * kFactors[7]);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i + 1), factor * m_h4s[i] * kFactors[8]);
        m_qp_problem.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * m_h4s[i] * kFactors[7]);
    }
    return;
}

auto RouteLineQuinticAccModel::setJerkCost(double jerk_weight) noexcept -> void {
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor =
            jerk_weight * m_reciprocal_h3s[i] * m_reciprocal_h2s[i] / m_time_scale;
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i), factor * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 1440.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i), factor * m_hs[i] * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * m_hs[i] * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i), factor * m_h2s[i] * 120.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * m_h2s[i] * 120.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * m_hs[i] * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * m_hs[i] * 720.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * m_h2s[i] * 120.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * m_h2s[i] * 120.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i), factor * m_h2s[i] * 192.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * m_h2s[i] * 336.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i), factor * m_h3s[i] * 72.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * m_h3s[i] * 48.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * m_h2s[i] * 192.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * m_h3s[i] * 48.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * m_h3s[i] * 72.0);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i), factor * m_h4s[i] * 9.0);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i + 1), -factor * m_h4s[i] * 6.0);
        m_qp_problem.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * m_h4s[i] * 9.0);
    }
    return;
}

auto RouteLineQuinticAccModel::setSnapCost(double snap_weight) noexcept -> void {
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor =
            snap_weight * m_reciprocal_h4s[i] * m_reciprocal_h3s[i] / m_time_scale;
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i), factor * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 86400.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i), factor * m_hs[i] * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * m_hs[i] * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i), factor * m_h2s[i] * 7200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * m_h2s[i] * 7200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * m_hs[i] * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * m_hs[i] * 43200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * m_h2s[i] * 7200.0);
        m_qp_problem.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * m_h2s[i] * 7200.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i), factor * m_h2s[i] * 10944.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * m_h2s[i] * 21312.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i), factor * m_h3s[i] * 3744.0);
        m_qp_problem.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * m_h3s[i] * 3456.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * m_h2s[i] * 10944.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * m_h3s[i] * 3456.0);
        m_qp_problem.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * m_h3s[i] * 3744.0);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i), factor * m_h4s[i] * 336.0);
        m_qp_problem.addQuadCostTerm(aIndex(i), aIndex(i + 1), -factor * m_h4s[i] * 528.0);
        m_qp_problem.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * m_h4s[i] * 336.0);
    }
    return;
}

auto RouteLineQuinticAccModel::solve() const
    -> std::pair<Motion1d, ::boyle::cvxopm::OsqpSolver::Info> {
    const ::boyle::cvxopm::OsqpSolver osqp_solver{m_settings};
    const auto [osqp_result, osqp_info] = osqp_solver.solve(m_qp_problem);
    const double v0{osqp_result.prim_vars[vIndex(0)]};
    const double vf{osqp_result.prim_vars[vIndex(m_num_samples - 1)]};
    const double a0{osqp_result.prim_vars[aIndex(0)]};
    const double af{osqp_result.prim_vars[aIndex(m_num_samples - 1)]};
    const std::array<Motion1d::BoundaryMode, 2> b0{
        Motion1d::BoundaryMode{1, v0}, Motion1d::BoundaryMode{2, a0}
    };
    const std::array<Motion1d::BoundaryMode, 2> bf{
        Motion1d::BoundaryMode{1, vf}, Motion1d::BoundaryMode{2, af}
    };
    std::vector<double> anchor_ss{
        osqp_result.prim_vars.cbegin(),
        osqp_result.prim_vars.cbegin() + sIndex(m_num_samples - 1) + 1
    };
    return std::make_pair(Motion1d{m_sample_ts, std::move(anchor_ss), b0, bf}, osqp_info);
}

auto RouteLineQuinticAccModel::clear() noexcept -> void {
    m_num_samples = 0;
    m_time_scale = 0.0;
    m_qp_problem.clear();
    m_sample_ts.clear();
    m_hs.clear();
    m_h2s.clear();
    m_h3s.clear();
    m_h4s.clear();
    m_reciprocal_hs.clear();
    m_reciprocal_h2s.clear();
    m_reciprocal_h3s.clear();
    m_reciprocal_h4s.clear();
    return;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
auto RouteLineQuinticAccModel::sIndex(int t_index) const noexcept -> int { return t_index; }

auto RouteLineQuinticAccModel::vIndex(int t_index) const noexcept -> int {
    return m_num_samples + t_index;
}

auto RouteLineQuinticAccModel::aIndex(int t_index) const noexcept -> int {
    return m_num_samples * 2 + t_index;
}

} // namespace boyle::kinetics
