/**
 * @file route_line_cubic_offset_model.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-01-10
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/kinetics/models/route_line_quintic_offset_model.hpp"

#include <array>
#include <format>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "boyle/common/utils/logging.hpp"
#include "boyle/math/functions/piecewise_functions/piecewise_linear_function1.hpp"
#include "boyle/math/quintic_interpolation.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::kinetics {

RouteLineQuinticOffsetModel::RouteLineQuinticOffsetModel(
    const std::vector<::boyle::math::Vec2d>& sketch_points, std::vector<double> sample_ss
) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
    if (sketch_points.size() < 2) {
        throw std::invalid_argument(std::format(
            "Invalid argument error detected: size of raw_sketch_points must larger than 2: "
            "sketch_points.size() = {0:d}.",
            sketch_points.size()
        ));
    }
    if (sample_ss.size() < 2) {
        throw std::invalid_argument(std::format(
            "Invalid argument error detected: size of sample_ss must larger than 2: "
            "sample_ss.size() = {0:d}.",
            sample_ss.size()
        ));
    }
#endif
    m_num_samples = static_cast<int>(sample_ss.size());
    m_length_scale = sample_ss.back() - sample_ss.front();
    m_qp_problem.resize(m_num_samples * 6, m_num_samples * 10);
    osqp_set_default_settings(&m_settings);
    m_settings.scaling = 0;
    m_sketch_curve = ::boyle::math::PiecewiseQuinticCurve2d{sketch_points};
    m_sample_points.reserve(m_num_samples);
    for (const double s : sample_ss) {
        m_sample_points.emplace_back(m_sketch_curve(s));
    }
    m_sample_ss = std::move(sample_ss);
    m_hs.reserve(m_num_samples - 1);
    m_h2s.reserve(m_num_samples - 1);
    m_h3s.reserve(m_num_samples - 1);
    m_h4s.reserve(m_num_samples - 1);
    m_h6s.reserve(m_num_samples - 1);
    m_h8s.reserve(m_num_samples - 1);
    m_reciprocal_hs.reserve(m_num_samples - 1);
    for (std::vector<double>::const_iterator it{m_sample_ss.cbegin() + 1}; it != m_sample_ss.cend();
         ++it) {
        const double diff = *it - *(it - 1);
        const double reciprocal_diff = 1.0 / diff;
        m_hs.push_back(diff);
        m_h2s.push_back(diff * diff);
        m_h3s.push_back(m_h2s.back() * diff);
        m_h4s.push_back(m_h2s.back() * m_h2s.back());
        m_h6s.push_back(m_h4s.back() * m_h2s.back());
        m_h8s.push_back(m_h4s.back() * m_h4s.back());
        m_reciprocal_hs.push_back(reciprocal_diff);
    }
    setIntegrationRelation();
}

RouteLineQuinticOffsetModel::RouteLineQuinticOffsetModel(
    const std::vector<::boyle::math::Vec2d>& raw_sketch_points, std::vector<double> sample_ss,
    const ::boyle::cvxopm::OsqpSolver::Settings& settings
) noexcept(!BOYLE_CHECK_PARAMS)
    : RouteLineQuinticOffsetModel{raw_sketch_points, std::move(sample_ss)} {
    m_settings = settings;
}

auto RouteLineQuinticOffsetModel::num_samples() const noexcept -> std::size_t {
    return m_sample_ss.size();
}

auto RouteLineQuinticOffsetModel::qp_problem() const noexcept
    -> const ::boyle::cvxopm::QpProblem<double, int>& {
    return m_qp_problem;
}

auto RouteLineQuinticOffsetModel::settings() const noexcept
    -> const ::boyle::cvxopm::OsqpSolver::Settings& {
    return m_settings;
}

auto RouteLineQuinticOffsetModel::setIntegrationRelation() noexcept -> void {
    for (int i{1}; i < m_num_samples - 1; ++i) {
        std::unordered_map<int, double> constrain_vec{
            {xIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
            {xIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
            {xIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
            {ddxIndex(i - 1), m_hs[i - 1]},
            {ddxIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
            {ddxIndex(i + 1), m_hs[i]},
            {d4xIndex(i - 1), -m_h3s[i - 1] * 7.0 / 60.0},
            {d4xIndex(i), -(m_h3s[i - 1] + m_h3s[i]) * 2.0 / 15.0},
            {d4xIndex(i + 1), -m_h3s[i] * 7.0 / 60.0}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 6 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
        constrain_vec = {
            {ddxIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
            {ddxIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
            {ddxIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
            {d4xIndex(i - 1), m_hs[i - 1]},
            {d4xIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
            {d4xIndex(i + 1), m_hs[i]}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 7 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
        constrain_vec = {
            {yIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
            {yIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
            {yIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
            {ddyIndex(i - 1), m_hs[i - 1]},
            {ddyIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
            {ddyIndex(i + 1), m_hs[i]},
            {d4yIndex(i - 1), -m_h3s[i - 1] * 7.0 / 60.0},
            {d4yIndex(i), -(m_h3s[i - 1] + m_h3s[i]) * 2.0 / 15.0},
            {d4yIndex(i + 1), -m_h3s[i] * 7.0 / 60.0}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 8 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
        constrain_vec = {
            {ddyIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
            {ddyIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
            {ddyIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
            {d4yIndex(i - 1), m_hs[i - 1]},
            {d4yIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
            {d4yIndex(i + 1), m_hs[i]}
        };
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 9 + i, constrain_vec, -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
    }
    return;
}

auto RouteLineQuinticOffsetModel::setHardBorders(const std::vector<HardBorder2d>& hard_borders
) noexcept -> void {
    std::vector<double> x_lower_bound(m_num_samples, std::numeric_limits<double>::lowest());
    std::vector<double> x_upper_bound(m_num_samples, std::numeric_limits<double>::max());
    std::vector<double> y_lower_bound(m_num_samples, std::numeric_limits<double>::lowest());
    std::vector<double> y_upper_bound(m_num_samples, std::numeric_limits<double>::max());
    for (const HardBorder2d& hard_border : hard_borders) {
        const int istart =
            ::boyle::math::nearestUpperElement(
                m_sample_points.cbegin(), m_sample_points.cend(), hard_border.bound_points.front()
            ) -
            m_sample_points.cbegin();
        const int iend =
            ::boyle::math::nearestUpperElement(
                m_sample_points.cbegin(), m_sample_points.cend(), hard_border.bound_points.back()
            ) -
            m_sample_points.cbegin();
        if (istart == m_num_samples) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! This soft border is not in the same region as "
                "the sketch points."
            );
            break;
        }
        if (iend == 0) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! This soft border is not in the same region as "
                "the sketch points."
            );
            break;
        }

        const std::vector<::boyle::math::Vec2d>& bound_points{hard_border.bound_points};
        std::vector<double> bound_ss;
        bound_ss.reserve(bound_points.size());
        for (const ::boyle::math::Vec2d& bound_point : bound_points) {
            bound_ss.push_back(m_sketch_curve.inverse(bound_point).s());
        }

        const ::boyle::math::PiecewiseLinearFunction1<::boyle::math::Vec2d, double> bound_func{
            bound_ss, bound_points
        };

        if (hard_border.chirality == ::boyle::kinetics::Chirality::RIGHT) {
            if (istart != 0) {
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                const double ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constraint_vec{
                    {xIndex(istart - 1), proration_coeffs[0]},
                    {xIndex(istart), proration_coeffs[1]},
                    {ddxIndex(istart - 1), proration_coeffs[2]},
                    {ddxIndex(istart), proration_coeffs[3]},
                    {d4xIndex(istart - 1), proration_coeffs[4]},
                    {d4xIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constraint_vec, bound_points.front().x(), std::numeric_limits<double>::max()
                );
                constraint_vec = {{yIndex(istart - 1), proration_coeffs[0]},
                                  {yIndex(istart), proration_coeffs[1]},
                                  {ddyIndex(istart - 1), proration_coeffs[2]},
                                  {ddyIndex(istart), proration_coeffs[3]},
                                  {d4xIndex(istart - 1), proration_coeffs[4]},
                                  {d4xIndex(istart), proration_coeffs[5]}};
                m_qp_problem.addConstrainTerm(
                    constraint_vec, bound_points.front().y(), std::numeric_limits<double>::max()
                );
            }
            for (int i{istart}; i < iend; ++i) {
                const ::boyle::math::Vec2d bound_point = bound_func(m_sample_ss[i]);
                x_lower_bound[i] = std::max(bound_point.x(), x_lower_bound[i]);
                y_lower_bound[i] = std::max(bound_point.y(), y_lower_bound[i]);
            }
            if (iend != m_num_samples) {
                const double h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                const double ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(iend - 1), proration_coeffs[0]},
                    {xIndex(iend), proration_coeffs[1]},
                    {ddxIndex(iend - 1), proration_coeffs[2]},
                    {ddxIndex(iend), proration_coeffs[3]},
                    {d4xIndex(iend - 1), proration_coeffs[4]},
                    {d4xIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, bound_points.back().x(), std::numeric_limits<double>::max()
                );
                constrain_vec = {{yIndex(iend - 1), proration_coeffs[0]},
                                 {yIndex(iend), proration_coeffs[1]},
                                 {ddyIndex(iend - 1), proration_coeffs[2]},
                                 {ddyIndex(iend), proration_coeffs[3]},
                                 {d4yIndex(iend - 1), proration_coeffs[4]},
                                 {d4yIndex(iend), proration_coeffs[5]}};
                m_qp_problem.addConstrainTerm(
                    constrain_vec, bound_points.back().y(), std::numeric_limits<double>::max()
                );
            }
        } else {
            if (istart != 0) {
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                const double ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(istart - 1), proration_coeffs[0]},
                    {xIndex(istart), proration_coeffs[1]},
                    {ddxIndex(istart - 1), proration_coeffs[2]},
                    {ddxIndex(istart), proration_coeffs[3]},
                    {d4xIndex(istart - 1), proration_coeffs[4]},
                    {d4xIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(), bound_points.front().x()
                );
                constrain_vec = {{yIndex(istart - 1), proration_coeffs[0]},
                                 {yIndex(istart), proration_coeffs[1]},
                                 {ddyIndex(istart - 1), proration_coeffs[2]},
                                 {ddyIndex(istart), proration_coeffs[3]},
                                 {d4yIndex(istart - 1), proration_coeffs[4]},
                                 {d4yIndex(istart), proration_coeffs[5]}};
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(), bound_points.front().y()
                );
            }
            for (int i{istart}; i < iend; ++i) {
                const ::boyle::math::Vec2d bound_point = bound_func(m_sample_ss[i]);
                x_upper_bound[i] = std::min(bound_point.x(), x_upper_bound[i]);
                y_upper_bound[i] = std::min(bound_point.y(), y_upper_bound[i]);
            }
            if (iend != m_num_samples) {
                const double h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                const double ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(iend - 1), proration_coeffs[0]},
                    {xIndex(iend), proration_coeffs[1]},
                    {ddxIndex(iend - 1), proration_coeffs[2]},
                    {ddxIndex(iend), proration_coeffs[3]},
                    {d4xIndex(iend - 1), proration_coeffs[4]},
                    {d4xIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(), bound_points.back().x()
                );
                constrain_vec = {{yIndex(iend - 1), proration_coeffs[0]},
                                 {yIndex(iend), proration_coeffs[1]},
                                 {ddyIndex(iend - 1), proration_coeffs[2]},
                                 {ddyIndex(iend), proration_coeffs[3]},
                                 {d4yIndex(iend - 1), proration_coeffs[4]},
                                 {d4yIndex(iend), proration_coeffs[5]}};
                m_qp_problem.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(), bound_points.back().y()
                );
            }
        }
    }

    for (int i{1}; i < m_num_samples - 1; ++i) {
        m_qp_problem.updateConstrainTerm(
            xIndex(i), {{xIndex(i), 1.0}}, x_lower_bound[i], x_upper_bound[i]
        );
        m_qp_problem.updateConstrainTerm(
            yIndex(i), {{yIndex(i), 1.0}}, y_lower_bound[i], y_upper_bound[i]
        );
    }
    return;
}

auto RouteLineQuinticOffsetModel::setSoftBorders(const std::vector<SoftBorder2d>& soft_borders
) noexcept -> void {
    for (const SoftBorder2d& soft_border : soft_borders) {
        const int istart =
            ::boyle::math::nearestUpperElement(
                m_sample_points.cbegin(), m_sample_points.cend(), soft_border.bound_points.front()
            ) -
            m_sample_points.cbegin();
        const int iend =
            ::boyle::math::nearestUpperElement(
                m_sample_points.cbegin(), m_sample_points.cend(), soft_border.bound_points.back()
            ) -
            m_sample_points.cbegin();
        if (istart == m_num_samples) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! This soft border is not in the same region as "
                "the sketch points."
            );
            break;
        }
        if (iend == 0) {
            BOYLE_LOG_WARN(
                "Invalid argument issue detected! This soft border is not in the same region as "
                "the sketch points."
            );
            break;
        }

        const std::vector<::boyle::math::Vec2d>& bound_points{soft_border.bound_points};
        std::vector<double> bound_ss;
        bound_ss.reserve(bound_points.size());
        for (const ::boyle::math::Vec2d& bound_point : bound_points) {
            bound_ss.push_back(m_sketch_curve.inverse(bound_point).s());
        }

        const ::boyle::math::PiecewiseLinearFunction1<::boyle::math::Vec2d, double> bound_func{
            bound_ss, bound_points
        };

        double factor{0.0};

        if (soft_border.chirality == ::boyle::kinetics::Chirality::RIGHT) {
            if (istart == iend) {
                factor = (bound_ss.back() - bound_ss.front()) / m_length_scale * 0.5;
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                double ratio{(m_sample_ss.front() - m_sample_ss[istart - 1]) / h};
                std::array<double, 6> proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(istart - 1), -proration_coeffs[0]},
                    {xIndex(istart), -proration_coeffs[1]},
                    {ddxIndex(istart - 1), -proration_coeffs[2]},
                    {ddxIndex(istart), -proration_coeffs[3]},
                    {d4xIndex(istart - 1), -proration_coeffs[4]},
                    {d4xIndex(istart), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.front().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(istart - 1), -proration_coeffs[0]},
                                 {yIndex(istart), -proration_coeffs[1]},
                                 {ddyIndex(istart - 1), -proration_coeffs[2]},
                                 {ddyIndex(istart), -proration_coeffs[3]},
                                 {d4yIndex(istart - 1), -proration_coeffs[4]},
                                 {d4yIndex(istart), -proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.front().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                ratio = (bound_ss.back() - m_sample_ss[iend - 1]) / h;
                proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                constrain_vec = {{xIndex(iend - 1), -proration_coeffs[0]},
                                 {xIndex(iend), -proration_coeffs[1]},
                                 {ddxIndex(iend - 1), -proration_coeffs[2]},
                                 {ddxIndex(iend), -proration_coeffs[3]},
                                 {d4xIndex(iend - 1), -proration_coeffs[4]},
                                 {d4xIndex(iend), -proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.back().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(iend - 1), -proration_coeffs[0]},
                                 {yIndex(iend), -proration_coeffs[1]},
                                 {ddyIndex(iend - 1), -proration_coeffs[2]},
                                 {ddyIndex(iend), -proration_coeffs[3]},
                                 {d4yIndex(iend - 1), -proration_coeffs[4]},
                                 {d4yIndex(iend), -proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.back().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                factor = (m_sample_ss[istart] - bound_ss.front()) / m_length_scale;
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                const double ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constraint_vec{
                    {xIndex(istart - 1), -proration_coeffs[0]},
                    {xIndex(istart), -proration_coeffs[1]},
                    {ddxIndex(istart - 1), -proration_coeffs[2]},
                    {ddxIndex(istart), -proration_coeffs[3]},
                    {d4xIndex(istart - 1), -proration_coeffs[4]},
                    {d4xIndex(istart), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constraint_vec), -bound_points.front().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constraint_vec = {{yIndex(istart - 1), -proration_coeffs[0]},
                                  {yIndex(istart), -proration_coeffs[1]},
                                  {ddyIndex(istart - 1), -proration_coeffs[2]},
                                  {ddyIndex(istart), -proration_coeffs[3]},
                                  {d4yIndex(istart - 1), -proration_coeffs[4]},
                                  {d4yIndex(istart), -proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constraint_vec), -bound_points.front().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
            }
            for (int i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = m_hs[i] / m_length_scale * 0.5;
                } else if (i == iend - 1) {
                    factor = m_hs[i - 1] / m_length_scale * 0.5;
                } else {
                    factor = (m_hs[i - 1] + m_hs[i]) / m_length_scale * 0.5;
                }
                const ::boyle::math::Vec2d bound_point = bound_func(m_sample_ss[i]);
                m_qp_problem.addClampCostTerm(
                    {{xIndex(i), -1.0}}, -bound_point.x(), soft_border.linear_weight * factor,
                    soft_border.quadratic_weight * factor
                );
                m_qp_problem.addClampCostTerm(
                    {{yIndex(i), -1.0}}, -bound_point.y(), soft_border.linear_weight * factor,
                    soft_border.quadratic_weight * factor
                );
            }
            if (iend != m_num_samples) {
                factor = (bound_ss.back() - m_sample_ss[iend - 1]) / m_length_scale;
                const double h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                const double ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(iend - 1), -proration_coeffs[0]},
                    {xIndex(iend), -proration_coeffs[1]},
                    {ddxIndex(iend - 1), -proration_coeffs[2]},
                    {ddxIndex(iend), -proration_coeffs[3]},
                    {d4xIndex(iend - 1), -proration_coeffs[4]},
                    {d4xIndex(iend), -proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.back().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(iend - 1), -proration_coeffs[0]},
                                 {yIndex(iend), -proration_coeffs[1]},
                                 {ddyIndex(iend - 1), -proration_coeffs[2]},
                                 {ddyIndex(iend), -proration_coeffs[3]},
                                 {d4yIndex(iend - 1), -proration_coeffs[4]},
                                 {d4yIndex(iend), -proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), -bound_points.back().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
            }
        } else {
            if (istart == iend) {
                factor = (bound_ss.back() - bound_ss.front()) / m_length_scale * 0.5;
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                double ratio{(m_sample_ss.front() - m_sample_ss[istart - 1]) / h};
                std::array<double, 6> proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(istart - 1), proration_coeffs[0]},
                    {xIndex(istart), proration_coeffs[1]},
                    {ddxIndex(istart - 1), proration_coeffs[2]},
                    {ddxIndex(istart), proration_coeffs[3]},
                    {d4xIndex(istart - 1), proration_coeffs[4]},
                    {d4xIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.front().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(istart - 1), proration_coeffs[0]},
                                 {yIndex(istart), proration_coeffs[1]},
                                 {ddyIndex(istart - 1), proration_coeffs[2]},
                                 {ddyIndex(istart), proration_coeffs[3]},
                                 {d4yIndex(istart - 1), proration_coeffs[4]},
                                 {d4yIndex(istart), proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.front().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                ratio = (bound_ss.back() - m_sample_ss[iend - 1]) / h;
                proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                constrain_vec = {{xIndex(iend - 1), proration_coeffs[0]},
                                 {xIndex(iend), proration_coeffs[1]},
                                 {ddxIndex(iend - 1), proration_coeffs[2]},
                                 {ddxIndex(iend), proration_coeffs[3]},
                                 {d4xIndex(iend - 1), proration_coeffs[4]},
                                 {d4xIndex(iend), proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.back().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(iend - 1), proration_coeffs[0]},
                                 {yIndex(iend), proration_coeffs[1]},
                                 {ddyIndex(iend - 1), proration_coeffs[2]},
                                 {ddyIndex(iend), proration_coeffs[3]},
                                 {d4yIndex(iend - 1), proration_coeffs[4]},
                                 {d4yIndex(iend), proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.back().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                factor = (m_sample_ss[istart] - bound_ss.front()) / m_length_scale;
                const double h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                const double ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(istart - 1), proration_coeffs[0]},
                    {xIndex(istart), proration_coeffs[1]},
                    {ddxIndex(istart - 1), proration_coeffs[2]},
                    {ddxIndex(istart), proration_coeffs[3]},
                    {d4xIndex(istart - 1), proration_coeffs[4]},
                    {d4xIndex(istart), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.front().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(istart - 1), proration_coeffs[0]},
                                 {yIndex(istart), proration_coeffs[1]},
                                 {ddyIndex(istart - 1), proration_coeffs[2]},
                                 {ddyIndex(istart), proration_coeffs[3]},
                                 {d4yIndex(istart - 1), proration_coeffs[4]},
                                 {d4yIndex(istart), proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.front().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
            }
            for (int i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = m_hs[i] / m_length_scale * 0.5;
                } else if (i == iend - 1) {
                    factor = m_hs[i - 1] / m_length_scale * 0.5;
                } else {
                    factor = (m_hs[i - 1] + m_hs[i]) / m_length_scale * 0.5;
                }
                const ::boyle::math::Vec2d bound_point = bound_func(m_sample_ss[i]);
                m_qp_problem.addClampCostTerm(
                    {{xIndex(i), 1.0}}, bound_point.x(), soft_border.linear_weight * factor,
                    soft_border.quadratic_weight * factor
                );
                m_qp_problem.addClampCostTerm(
                    {{yIndex(i), 1.0}}, bound_point.y(), soft_border.linear_weight * factor,
                    soft_border.quadratic_weight * factor
                );
            }
            if (iend != m_num_samples) {
                factor = (bound_ss.back() - m_sample_ss[iend - 1]) / m_length_scale;
                const double h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                const double ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs =
                    ::boyle::math::quinerpCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {xIndex(iend - 1), proration_coeffs[0]},
                    {xIndex(iend), proration_coeffs[1]},
                    {ddxIndex(iend - 1), proration_coeffs[2]},
                    {ddxIndex(iend), proration_coeffs[3]},
                    {d4xIndex(iend - 1), proration_coeffs[4]},
                    {d4xIndex(iend), proration_coeffs[5]}
                };
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.back().x(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
                constrain_vec = {{yIndex(iend - 1), proration_coeffs[0]},
                                 {yIndex(iend), proration_coeffs[1]},
                                 {ddyIndex(iend - 1), proration_coeffs[2]},
                                 {ddyIndex(iend), proration_coeffs[3]},
                                 {d4yIndex(iend - 1), proration_coeffs[4]},
                                 {d4yIndex(iend), proration_coeffs[5]}};
                m_qp_problem.addClampCostTerm(
                    std::move(constrain_vec), bound_points.back().y(),
                    soft_border.linear_weight * factor, soft_border.quadratic_weight * factor
                );
            }
        }
    }
    return;
}

auto RouteLineQuinticOffsetModel::setDdxRange(double ddx_min, double ddx_max) noexcept -> void {
    for (int i{0}; i < m_num_samples; i++) {
        m_qp_problem.updateConstrainTerm(ddxIndex(i), {{ddxIndex(i), 1.0}}, ddx_min, ddx_max);
    }
    return;
}

auto RouteLineQuinticOffsetModel::setDdyRange(double ddy_min, double ddy_max) noexcept -> void {
    for (int i{0}; i < m_num_samples; i++) {
        m_qp_problem.updateConstrainTerm(ddyIndex(i), {{ddyIndex(i), 1.0}}, ddy_min, ddy_max);
    }
    return;
}

auto RouteLineQuinticOffsetModel::setD4xRange(double ddx_min, double ddx_max) noexcept -> void {
    for (int i{0}; i < m_num_samples; i++) {
        m_qp_problem.updateConstrainTerm(d4xIndex(i), {{d4xIndex(i), 1.0}}, ddx_min, ddx_max);
    }
    return;
}

auto RouteLineQuinticOffsetModel::setD4yRange(double ddy_min, double ddy_max) noexcept -> void {
    for (int i{0}; i < m_num_samples; i++) {
        m_qp_problem.updateConstrainTerm(d4yIndex(i), {{d4yIndex(i), 1.0}}, ddy_min, ddy_max);
    }
    return;
}

auto RouteLineQuinticOffsetModel::setInitialState(
    ::boyle::math::Vec2d r0, ::boyle::math::Vec2d t0, ::boyle::math::Vec2d n0,
    ::boyle::math::Vec2d j0
) noexcept -> void {
    m_qp_problem.updateConstrainTerm(
        xIndex(0), {{xIndex(0), 1.0}}, r0.x() - ::boyle::math::kEpsilon,
        r0.x() + ::boyle::math::kEpsilon
    );
    if (!std::isnan(t0.x())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 6,
            {{xIndex(0), -m_reciprocal_hs[0]},
             {xIndex(1), m_reciprocal_hs[0]},
             {ddxIndex(0), -m_hs[0] / 3.0},
             {ddxIndex(1), -m_hs[0] / 6.0},
             {d4xIndex(0), m_h3s[0] / 45.0},
             {d4xIndex(1), m_h3s[0] * 7.0 / 360.0}},
            t0.x() - ::boyle::math::kEpsilon, t0.x() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(n0.x())) {
        m_qp_problem.updateConstrainTerm(
            ddxIndex(0), {{ddxIndex(0), 1.0}}, n0.x() - ::boyle::math::kEpsilon,
            n0.x() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(j0.x())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 7,
            {{ddxIndex(0), -m_reciprocal_hs[0]},
             {ddxIndex(1), m_reciprocal_hs[0]},
             {d4xIndex(0), -m_hs[0] / 3.0},
             {d4xIndex(1), -m_hs[0] / 6.0}},
            j0.x() - ::boyle::math::kEpsilon, j0.x() + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        d4xIndex(0), {{d4xIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::max()
    );*/

    m_qp_problem.updateConstrainTerm(
        yIndex(0), {{yIndex(0), 1.0}}, r0.y() - ::boyle::math::kEpsilon,
        r0.y() + ::boyle::math::kEpsilon
    );
    if (!std::isnan(t0.y())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 8,
            {{yIndex(0), -m_reciprocal_hs[0]},
             {yIndex(1), m_reciprocal_hs[0]},
             {ddyIndex(0), -m_hs[0] / 3.0},
             {ddyIndex(1), -m_hs[0] / 6.0},
             {d4yIndex(0), m_h3s[0] / 45.0},
             {d4yIndex(1), m_h3s[0] * 7.0 / 360.0}},
            t0.y() - ::boyle::math::kEpsilon, t0.y() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(n0.y())) {
        m_qp_problem.updateConstrainTerm(
            ddyIndex(0), {{ddyIndex(0), 1.0}}, n0.y() - ::boyle::math::kEpsilon,
            n0.y() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(j0.y())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 9,
            {{ddyIndex(0), -m_reciprocal_hs[0]},
             {ddyIndex(1), m_reciprocal_hs[0]},
             {d4yIndex(0), -m_hs[0] / 3.0},
             {d4yIndex(1), -m_hs[0] / 6.0}},
            j0.y() - ::boyle::math::kEpsilon, j0.y() + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        d4yIndex(0), {{d4yIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::max()
    );*/
    return;
}

auto RouteLineQuinticOffsetModel::setFinalState(
    ::boyle::math::Vec2d rf, ::boyle::math::Vec2d tf, ::boyle::math::Vec2d nf,
    ::boyle::math::Vec2d jf
) noexcept -> void {
    m_qp_problem.updateConstrainTerm(
        xIndex(m_num_samples - 1), {{xIndex(m_num_samples - 1), 1.0}},
        rf.x() - ::boyle::math::kEpsilon, rf.x() + ::boyle::math::kEpsilon
    );
    if (!std::isnan(tf.x())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 7 - 1,
            {{xIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
             {xIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
             {ddxIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
             {ddxIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0},
             {d4xIndex(m_num_samples - 2), -m_h3s[m_num_samples - 2] * 7.0 / 360.0},
             {d4xIndex(m_num_samples - 1), -m_h3s[m_num_samples - 2] / 45.0}},
            tf.x() - ::boyle::math::kEpsilon, tf.x() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(nf.x())) {
        m_qp_problem.updateConstrainTerm(
            ddxIndex(m_num_samples - 1), {{ddxIndex(m_num_samples - 1), 1.0}},
            nf.x() - ::boyle::math::kEpsilon, nf.x() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(jf.x())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 8 - 1,
            {{ddxIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
             {ddxIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
             {d4xIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
             {d4xIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0}},
            jf.x() - ::boyle::math::kEpsilon, jf.x() + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        d4xIndex(m_num_samples - 1), {{d4xIndex(m_num_samples - 1), 1.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );*/

    m_qp_problem.updateConstrainTerm(
        yIndex(m_num_samples - 1), {{yIndex(m_num_samples - 1), 1.0}},
        rf.y() - ::boyle::math::kEpsilon, rf.y() + ::boyle::math::kEpsilon
    );
    if (!std::isnan(tf.y())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 9 - 1,
            {{yIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
             {yIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
             {ddyIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
             {ddyIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0},
             {d4yIndex(m_num_samples - 2), -m_h3s[m_num_samples - 2] * 7.0 / 360.0},
             {d4yIndex(m_num_samples - 1), -m_h3s[m_num_samples - 2] / 45.0}},
            tf.y() - ::boyle::math::kEpsilon, tf.y() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(nf.y())) {
        m_qp_problem.updateConstrainTerm(
            ddyIndex(m_num_samples - 1), {{ddyIndex(m_num_samples - 1), 1.0}},
            nf.y() - ::boyle::math::kEpsilon, nf.y() + ::boyle::math::kEpsilon
        );
    }
    if (!std::isnan(jf.y())) {
        m_qp_problem.updateConstrainTerm(
            m_num_samples * 8 - 1,
            {{ddyIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
             {ddyIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
             {d4yIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
             {d4yIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0}},
            jf.y() - ::boyle::math::kEpsilon, jf.y() + ::boyle::math::kEpsilon
        );
    }
    /*m_qp_problem.updateConstrainTerm(
        d4yIndex(m_num_samples - 1), {{d4yIndex(m_num_samples - 1), 1.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );*/
    return;
}

auto RouteLineQuinticOffsetModel::setOffsetCost(double offset_weight) noexcept -> void {
    constexpr std::array<double, 10> kFactors{1.0 / 3.0,       2.0 / 45.0,       7.0 / 180.0,
                                              2.0 / 945.0,     4.0 / 945.0,      31.0 / 7560.0,
                                              2.0 / 4725.0,    127.0 / 302400.0, 2.0 / 93555.0,
                                              73.0 / 1710720.0};
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor = offset_weight * m_hs[i] / m_length_scale;

        m_qp_problem.addQuadCostTerm(xIndex(i), xIndex(i), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(xIndex(i), xIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(xIndex(i), ddxIndex(i), -factor * m_h2s[i] * kFactors[1]);
        m_qp_problem.addQuadCostTerm(xIndex(i), ddxIndex(i + 1), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(xIndex(i), d4xIndex(i), factor * m_h4s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(xIndex(i), d4xIndex(i + 1), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(xIndex(i + 1), xIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(xIndex(i + 1), ddxIndex(i), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(
            xIndex(i + 1), ddxIndex(i + 1), -factor * m_h2s[i] * kFactors[1]
        );
        m_qp_problem.addQuadCostTerm(xIndex(i + 1), d4xIndex(i), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(
            xIndex(i + 1), d4xIndex(i + 1), factor * m_h4s[i] * kFactors[4]
        );
        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(ddxIndex(i), d4xIndex(i), -factor * m_h6s[i] * kFactors[6]);
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i), d4xIndex(i + 1), -factor * m_h6s[i] * kFactors[7]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[3]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), d4xIndex(i), -factor * m_h6s[i] * kFactors[7]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), d4xIndex(i + 1), -factor * m_h6s[i] * kFactors[6]
        );
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i), factor * m_h8s[i] * kFactors[8]);
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[9]);
        m_qp_problem.addQuadCostTerm(
            d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[8]
        );
        m_qp_problem.addLinCostTerm(
            xIndex(i),
            -factor * (m_sample_points[i].x() * 2.0 + m_sample_points[i + 1].x()) * kFactors[0]
        );
        m_qp_problem.addLinCostTerm(
            xIndex(i + 1),
            -factor * (m_sample_points[i].x() + m_sample_points[i + 1].x() * 2.0) * kFactors[0]
        );
        m_qp_problem.addLinCostTerm(
            ddxIndex(i),
            factor * m_h2s[i] *
                (m_sample_points[i].x() * kFactors[1] + m_sample_points[i + 1].x() * kFactors[2])
        );
        m_qp_problem.addLinCostTerm(
            ddxIndex(i + 1),
            factor * m_h2s[i] *
                (m_sample_points[i].x() * kFactors[2] + m_sample_points[i + 1].x() * kFactors[1])
        );
        m_qp_problem.addLinCostTerm(
            d4xIndex(i),
            -factor * m_h4s[i] *
                (m_sample_points[i].x() * kFactors[4] + m_sample_points[i + 1].x() * kFactors[5])
        );
        m_qp_problem.addLinCostTerm(
            d4xIndex(i + 1),
            -factor * m_h4s[i] *
                (m_sample_points[i].x() * kFactors[5] + m_sample_points[i + 1].x() * kFactors[4])
        );

        m_qp_problem.addQuadCostTerm(yIndex(i), yIndex(i), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(yIndex(i), yIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(yIndex(i), ddyIndex(i), -factor * m_h2s[i] * kFactors[1]);
        m_qp_problem.addQuadCostTerm(yIndex(i), ddyIndex(i + 1), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(yIndex(i), d4yIndex(i), factor * m_h4s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(yIndex(i), d4yIndex(i + 1), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(yIndex(i + 1), yIndex(i + 1), factor * kFactors[0]);
        m_qp_problem.addQuadCostTerm(yIndex(i + 1), ddyIndex(i), -factor * m_h2s[i] * kFactors[2]);
        m_qp_problem.addQuadCostTerm(
            yIndex(i + 1), ddyIndex(i + 1), -factor * m_h2s[i] * kFactors[1]
        );
        m_qp_problem.addQuadCostTerm(yIndex(i + 1), d4yIndex(i), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(
            yIndex(i + 1), d4yIndex(i + 1), factor * m_h4s[i] * kFactors[4]
        );
        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[5]);
        m_qp_problem.addQuadCostTerm(ddyIndex(i), d4yIndex(i), -factor * m_h6s[i] * kFactors[6]);
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i), d4yIndex(i + 1), -factor * m_h6s[i] * kFactors[7]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[3]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), d4yIndex(i), -factor * m_h6s[i] * kFactors[7]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), d4yIndex(i + 1), -factor * m_h6s[i] * kFactors[6]
        );
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i), factor * m_h8s[i] * kFactors[8]);
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[9]);
        m_qp_problem.addQuadCostTerm(
            d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[8]
        );
        m_qp_problem.addLinCostTerm(
            yIndex(i),
            -factor * (m_sample_points[i].y() * 2.0 + m_sample_points[i + 1].y()) * kFactors[0]
        );
        m_qp_problem.addLinCostTerm(
            yIndex(i + 1),
            -factor * (m_sample_points[i].y() + m_sample_points[i + 1].y() * 2.0) * kFactors[0]
        );
        m_qp_problem.addLinCostTerm(
            ddyIndex(i),
            factor * m_h2s[i] *
                (m_sample_points[i].y() * kFactors[1] + m_sample_points[i + 1].y() * kFactors[2])
        );
        m_qp_problem.addLinCostTerm(
            ddyIndex(i + 1),
            factor * m_h2s[i] *
                (m_sample_points[i].y() * kFactors[2] + m_sample_points[i + 1].y() * kFactors[1])
        );
        m_qp_problem.addLinCostTerm(
            d4yIndex(i),
            -factor * m_h4s[i] *
                (m_sample_points[i].y() * kFactors[4] + m_sample_points[i + 1].y() * kFactors[5])
        );
        m_qp_problem.addLinCostTerm(
            d4yIndex(i + 1),
            -factor * m_h4s[i] *
                (m_sample_points[i].y() * kFactors[5] + m_sample_points[i + 1].y() * kFactors[4])
        );
    }
    return;
}

auto RouteLineQuinticOffsetModel::setCurvatureCost(double curvature_weight) noexcept -> void {
    constexpr std::array<double, 5> kFactors{
        1.0 / 3.0, 2.0 / 45.0, 7.0 / 180.0, 2.0 / 945.0, 31.0 / 7560.0
    };
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor = curvature_weight * m_reciprocal_hs[i] * m_reciprocal_hs[i] *
                              m_reciprocal_hs[i] / m_length_scale;

        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(ddxIndex(i), d4xIndex(i), -factor * m_h6s[i] * kFactors[1]);
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i), d4xIndex(i + 1), -factor * m_h6s[i] * kFactors[2]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[0]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), d4xIndex(i), -factor * m_h6s[i] * kFactors[2]
        );
        m_qp_problem.addQuadCostTerm(
            ddxIndex(i + 1), d4xIndex(i + 1), -factor * m_h6s[i] * kFactors[1]
        );
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i), factor * m_h8s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(
            d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[3]
        );

        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[0]);
        m_qp_problem.addQuadCostTerm(ddyIndex(i), d4yIndex(i), -factor * m_h6s[i] * kFactors[1]);
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i), d4yIndex(i + 1), -factor * m_h6s[i] * kFactors[2]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[0]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), d4yIndex(i), -factor * m_h6s[i] * kFactors[2]
        );
        m_qp_problem.addQuadCostTerm(
            ddyIndex(i + 1), d4yIndex(i + 1), -factor * m_h6s[i] * kFactors[1]
        );
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i), factor * m_h8s[i] * kFactors[3]);
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[4]);
        m_qp_problem.addQuadCostTerm(
            d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[3]
        );
    }
    return;
}

auto RouteLineQuinticOffsetModel::setDCurvatureCost(double dcurvature_weight) noexcept -> void {
    for (int i{0}; i < m_num_samples - 1; ++i) {
        const double factor = dcurvature_weight * m_reciprocal_hs[i] * m_reciprocal_hs[i] *
                              m_reciprocal_hs[i] * m_reciprocal_hs[i] * m_reciprocal_hs[i] /
                              m_length_scale;
        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i]);
        m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i + 1), -factor * m_h4s[i] * 2.0);
        m_qp_problem.addQuadCostTerm(ddxIndex(i + 1), ddxIndex(i + 1), factor * m_h4s[i]);
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i), factor * m_h8s[i] / 45.0);
        m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * 7.0 / 180.0);
        m_qp_problem.addQuadCostTerm(d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] / 45.0);

        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i]);
        m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i + 1), -factor * m_h4s[i] * 2.0);
        m_qp_problem.addQuadCostTerm(ddyIndex(i + 1), ddyIndex(i + 1), factor * m_h4s[i]);
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i), factor * m_h8s[i] / 45.0);
        m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * 7.0 / 180.0);
        m_qp_problem.addQuadCostTerm(d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] / 45.0);
    }
    return;
}

auto RouteLineQuinticOffsetModel::solve() const
    -> std::pair<Path2d, ::boyle::cvxopm::OsqpSolver::Info> {
    const ::boyle::cvxopm::OsqpSolver osqp_solver{m_settings};
    const auto [osqp_result, osqp_info] = osqp_solver.solve(m_qp_problem);
    const double ddx0{osqp_result.prim_vars[ddxIndex(0)]};
    const double ddxf{osqp_result.prim_vars[ddxIndex(m_num_samples - 1)]};
    const double ddy0{osqp_result.prim_vars[ddyIndex(0)]};
    const double ddyf{osqp_result.prim_vars[ddyIndex(m_num_samples - 1)]};
    const double d4x0{osqp_result.prim_vars[d4xIndex(0)]};
    const double d4xf{osqp_result.prim_vars[d4xIndex(m_num_samples - 1)]};
    const double d4y0{osqp_result.prim_vars[d4yIndex(0)]};
    const double d4yf{osqp_result.prim_vars[d4yIndex(m_num_samples - 1)]};
    const std::array<Path2d::BoundaryMode, 2> b0{
        Path2d::BoundaryMode{2, {ddx0, ddy0}}, Path2d::BoundaryMode{4, {d4x0, d4y0}}
    };
    const std::array<Path2d::BoundaryMode, 2> bf{
        Path2d::BoundaryMode{2, {ddxf, ddyf}}, Path2d::BoundaryMode{4, {d4xf, d4yf}}
    };
    std::vector<::boyle::math::Vec2d> anchor_points = ::boyle::math::squeeze<::boyle::math::Vec2d>(
        {osqp_result.prim_vars.cbegin() + xIndex(0),
         osqp_result.prim_vars.cbegin() + xIndex(m_num_samples - 1) + 1},
        {osqp_result.prim_vars.cbegin() + yIndex(0),
         osqp_result.prim_vars.cbegin() + yIndex(m_num_samples - 1) + 1}
    );
    return std::make_pair(Path2d{std::move(anchor_points), b0, bf}, osqp_info);
}

auto RouteLineQuinticOffsetModel::clear() noexcept -> void {
    m_sketch_curve = ::boyle::math::PiecewiseQuinticCurve2d{};
    m_num_samples = 0;
    m_sample_points.clear();
    m_length_scale = 0.0;
    m_qp_problem.clear();
    m_sample_ss.clear();
    m_hs.clear();
    m_h2s.clear();
    m_h3s.clear();
    m_h4s.clear();
    m_h6s.clear();
    m_h8s.clear();
    m_reciprocal_hs.clear();
    return;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
auto RouteLineQuinticOffsetModel::xIndex(int s_index) const noexcept -> int { return s_index; }

auto RouteLineQuinticOffsetModel::ddxIndex(int s_index) const noexcept -> int {
    return m_num_samples + s_index;
}

auto RouteLineQuinticOffsetModel::d4xIndex(int s_index) const noexcept -> int {
    return m_num_samples * 2 + s_index;
}

auto RouteLineQuinticOffsetModel::yIndex(int s_index) const noexcept -> int {
    return m_num_samples * 3 + s_index;
}

auto RouteLineQuinticOffsetModel::ddyIndex(int s_index) const noexcept -> int {
    return m_num_samples * 4 + s_index;
}

auto RouteLineQuinticOffsetModel::d4yIndex(int s_index) const noexcept -> int {
    return m_num_samples * 5 + s_index;
}

} // namespace boyle::kinetics
