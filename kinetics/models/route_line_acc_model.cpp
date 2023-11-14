/**
 * @file route_line_acc_model.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "kinetics/models/route_line_acc_model.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "spdlog/spdlog.h"

namespace tiny_pnc {
namespace kinetics {

std::array<double, 6> RouteLineAccModel::prorationCoeffs(double ratio, double scale) {
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
        (ratio3 - ratio4 * 2.0 + ratio5) * scale2 * 0.5};
}

RouteLineAccModel::RouteLineAccModel(std::vector<double> sample_ts) {
    if (sample_ts.size() < 2) {
        std::string error_msg = std::format(
            "Invalid argument error detected: size of sample_ts must larger than 2: "
            "sample_ts.size() = {0:d}.",
            sample_ts.size()
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    num_samples_ = sample_ts.size();
    qp_problem_.resize(num_samples_ * 3, num_samples_ * 5);
    osqp_set_default_settings(&settings_);
    settings_.scaling = 0;
    sample_ts_ = std::move(sample_ts);
    time_scale_ = sample_ts_.back() - sample_ts_.front();
    hs_.reserve(num_samples_ - 1);
    h2s_.reserve(num_samples_ - 1);
    h3s_.reserve(num_samples_ - 1);
    h4s_.reserve(num_samples_ - 1);
    reciprocal_hs_.reserve(num_samples_ - 1);
    reciprocal_h2s_.reserve(num_samples_ - 1);
    reciprocal_h3s_.reserve(num_samples_ - 1);
    reciprocal_h4s_.reserve(num_samples_ - 1);
    for (std::vector<double>::const_iterator it = sample_ts_.cbegin() + 1; it != sample_ts_.cend();
         ++it) {
        const double diff = *it - *(it - 1);
        const double reciprocal_diff = 1.0 / diff;
        hs_.push_back(diff);
        h2s_.push_back(diff * diff);
        h3s_.push_back(h2s_.back() * diff);
        h4s_.push_back(h3s_.back() * diff);
        reciprocal_hs_.push_back(reciprocal_diff);
        reciprocal_h2s_.push_back(reciprocal_diff * reciprocal_diff);
        reciprocal_h3s_.push_back(reciprocal_h2s_.back() * reciprocal_diff);
        reciprocal_h4s_.push_back(reciprocal_h3s_.back() * reciprocal_diff);
    }
    setIntegrationRelation();
    setInitialState();
    setFinalState();
    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(
            sIndex(i), {{sIndex(i), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
        qp_problem_.updateConstrainTerm(
            vIndex(i), {{vIndex(i), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
        qp_problem_.updateConstrainTerm(
            aIndex(i), {{aIndex(i), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
}

RouteLineAccModel::RouteLineAccModel(
    std::vector<double> sample_ts, tiny_pnc::math::OsqpSettings settings
)
    : RouteLineAccModel(std::move(sample_ts)) {
    settings_ = settings;
}

const std::size_t& RouteLineAccModel::num_samples() const noexcept { return num_samples_; }

void RouteLineAccModel::setIntegrationRelation() noexcept {
    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        const std::unordered_map<int, double> constrain_vec_1{
            {sIndex(i - 1), reciprocal_h3s_[i - 1] * 20.0},
            {sIndex(i), -(reciprocal_h3s_[i - 1] + reciprocal_h3s_[i]) * 20.0},
            {sIndex(i + 1), reciprocal_h3s_[i] * 20.0},
            {vIndex(i - 1), reciprocal_h2s_[i - 1] * 8.0},
            {vIndex(i), (reciprocal_h2s_[i - 1] - reciprocal_h2s_[i]) * 12.0},
            {vIndex(i + 1), -reciprocal_h2s_[i] * 8.0},
            {aIndex(i - 1), reciprocal_hs_[i - 1]},
            {aIndex(i), -(reciprocal_hs_[i - 1] + reciprocal_hs_[i]) * 3.0},
            {aIndex(i + 1), reciprocal_hs_[i]}};
        const std::unordered_map<int, double> constrain_vec_2{
            {sIndex(i - 1), reciprocal_h4s_[i - 1] * 30.0},
            {sIndex(i), -(reciprocal_h4s_[i - 1] - reciprocal_h4s_[i]) * 30.0},
            {sIndex(i + 1), -reciprocal_h4s_[i] * 30.0},
            {vIndex(i - 1), reciprocal_h3s_[i - 1] * 14.0},
            {vIndex(i), (reciprocal_h3s_[i - 1] + reciprocal_h3s_[i]) * 16.0},
            {vIndex(i + 1), reciprocal_h3s_[i] * 14.0},
            {aIndex(i - 1), reciprocal_h2s_[i - 1] * 2.0},
            {aIndex(i), -(reciprocal_h2s_[i - 1] - reciprocal_h2s_[i]) * 3.0},
            {aIndex(i + 1), -reciprocal_h2s_[i] * 2.0}};
        qp_problem_.updateConstrainTerm(
            num_samples_ * 3 + i, constrain_vec_1, -tiny_pnc::math::kEpsilon,
            tiny_pnc::math::kEpsilon
        );
        qp_problem_.updateConstrainTerm(
            num_samples_ * 4 + i, constrain_vec_2, -tiny_pnc::math::kEpsilon,
            tiny_pnc::math::kEpsilon
        );
    }
    return;
}

void RouteLineAccModel::setSoftFences(const std::vector<SoftFence1d>& soft_fences) noexcept {
    for (const SoftFence1d& soft_fence : soft_fences) {
        const std::size_t istart =
            std::upper_bound(sample_ts_.cbegin(), sample_ts_.cend(), soft_fence.bound_line.minT()) -
            sample_ts_.cbegin();
        const std::size_t iend =
            std::upper_bound(sample_ts_.cbegin(), sample_ts_.cend(), soft_fence.bound_line.maxT()) -
            sample_ts_.cbegin();
        if (istart == num_samples_) {
            spdlog::warn(
                "Invalid argument issue detected! The minT() of hard_fence.bound_line should be "
                "less than sample_ts_.back(): hard_fence.bound_line.minT() = {0:f} while "
                "sample_ts_.back() = {1:f}.",
                soft_fence.bound_line.minT(), sample_ts_.back()
            );
            break;
        }
        if (iend == 0) {
            spdlog::warn(
                "Invalid argument issue detected! The maxT() of hard_fence.bound_line should be "
                "larger than sample_ts_.front(): hard_fence.bound_line.maxT() = {0:f} while "
                "sample_ts_.front() = {1:f}.",
                soft_fence.bound_line.maxT(), sample_ts_.front()
            );
            break;
        }

        double factor;

        if (soft_fence.actio == tiny_pnc::common::Actio::PUSHING) {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), -proration_coeffs[0]},
                    {sIndex(istart), -proration_coeffs[1]},
                    {vIndex(istart - 1), -proration_coeffs[2]},
                    {vIndex(istart), -proration_coeffs[3]},
                    {aIndex(istart - 1), -proration_coeffs[4]},
                    {aIndex(istart), -proration_coeffs[5]}};
                factor = (sample_ts_[istart] - soft_fence.bound_line.minT()) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (std::size_t i = istart; i < iend; ++i) {
                if (i == istart) {
                    factor = hs_[i] / time_scale_ * 0.5;
                } else if (i == iend - 1) {
                    factor = hs_[i - 1] / time_scale_ * 0.5;
                } else {
                    factor = (hs_[i - 1] + hs_[i]) / time_scale_ * 0.5;
                }
                qp_problem_.addClampCostTerm(
                    {{sIndex(i), -1.0}}, -soft_fence.bound_line(sample_ts_[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), -proration_coeffs[0]}, {sIndex(iend), -proration_coeffs[1]},
                    {vIndex(iend - 1), -proration_coeffs[2]}, {vIndex(iend), -proration_coeffs[3]},
                    {aIndex(iend - 1), -proration_coeffs[4]}, {aIndex(iend), -proration_coeffs[5]}};
                factor = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        } else {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}};
                factor = (sample_ts_[istart] - soft_fence.bound_line.minT()) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (std::size_t i = istart; i < iend; ++i) {
                if (i == istart) {
                    factor = hs_[i] / time_scale_ * 0.5;
                } else if (i == iend - 1) {
                    factor = hs_[i - 1] / time_scale_ * 0.5;
                } else {
                    factor = (hs_[i - 1] + hs_[i]) / time_scale_ * 0.5;
                }
                qp_problem_.addClampCostTerm(
                    {{sIndex(i), 1.0}}, soft_fence.bound_line(sample_ts_[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}};
                factor = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        }
    }
    return;
}

void RouteLineAccModel::setHardFences(const std::vector<HardFence1d>& hard_fences) noexcept {
    std::vector<double> lower_bound(num_samples_, std::numeric_limits<double>::lowest());
    std::vector<double> upper_bound(num_samples_, std::numeric_limits<double>::max());
    for (const HardFence1d& hard_fence : hard_fences) {
        const std::size_t istart =
            std::upper_bound(sample_ts_.cbegin(), sample_ts_.cend(), hard_fence.bound_line.minT()) -
            sample_ts_.cbegin();
        const std::size_t iend =
            std::upper_bound(sample_ts_.cbegin(), sample_ts_.cend(), hard_fence.bound_line.maxT()) -
            sample_ts_.cbegin();
        if (istart == num_samples_) {
            spdlog::warn(
                "Invalid argument issue detected! The minT() of hard_fence.bound_line should be "
                "less than sample_ts_.back(): hard_fence.bound_line.minT() = {0:f} while "
                "sample_ts_.back() = {1:f}.",
                hard_fence.bound_line.minT(), sample_ts_.back()
            );
            break;
        }
        if (iend == 0) {
            spdlog::warn(
                "Invalid argument issue detected! The maxT() of hard_fence.bound_line should be "
                "larger than sample_ts_.front(): hard_fence.bound_line.maxT() = {0:f} while "
                "sample_ts_.front() = {1:f}.",
                hard_fence.bound_line.maxT(), sample_ts_.front()
            );
            break;
        }

        if (hard_fence.actio == tiny_pnc::common::Actio::PUSHING) {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(hard_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, hard_fence.bound_line.ys().front(),
                    std::numeric_limits<double>::max()
                );
            }
            for (std::size_t i = istart; i < iend; ++i) {
                lower_bound[i] = std::max(hard_fence.bound_line(sample_ts_[i]), lower_bound[i]);
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(hard_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, hard_fence.bound_line.ys().back(),
                    std::numeric_limits<double>::max()
                );
            }
        } else {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(hard_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]},
                    {aIndex(istart - 1), proration_coeffs[4]},
                    {aIndex(istart), proration_coeffs[5]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(),
                    hard_fence.bound_line.ys().front()
                );
            }
            for (std::size_t i = istart; i < iend; ++i) {
                upper_bound[i] = std::min(hard_fence.bound_line(sample_ts_[i]), upper_bound[i]);
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(hard_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 6> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]}, {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]}, {vIndex(iend), proration_coeffs[3]},
                    {aIndex(iend - 1), proration_coeffs[4]}, {aIndex(iend), proration_coeffs[5]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(),
                    hard_fence.bound_line.ys().back()
                );
            }
        }
    }

    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(
            sIndex(i), {{sIndex(i), 1.0}}, lower_bound[i], upper_bound[i]
        );
    }
    return;
}

void RouteLineAccModel::setVelocityRange(double lower_bound, double upper_bound) noexcept {
    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(vIndex(i), {{vIndex(i), 1.0}}, lower_bound, upper_bound);
    }
    return;
}

void RouteLineAccModel::setAccelRange(double lower_bound, double upper_bound) noexcept {
    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(aIndex(i), {{aIndex(i), 1.0}}, lower_bound, upper_bound);
    }
    return;
}

void RouteLineAccModel::setInitialState(double s0, double v0, double a0) noexcept {
    if (!std::isnan(s0)) {
        qp_problem_.updateConstrainTerm(
            sIndex(0), {{sIndex(0), 1.0}}, s0 - tiny_pnc::math::kEpsilon,
            s0 + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            sIndex(0), {{sIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(v0)) {
        qp_problem_.updateConstrainTerm(
            vIndex(0), {{vIndex(0), 1.0}}, v0 - tiny_pnc::math::kEpsilon,
            v0 + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            vIndex(0), {{vIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(a0)) {
        qp_problem_.updateConstrainTerm(
            aIndex(0), {{aIndex(0), 1.0}}, a0 - tiny_pnc::math::kEpsilon,
            a0 + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            aIndex(0), {{aIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
    qp_problem_.updateConstrainTerm(
        num_samples_ * 3,
        {{sIndex(0), -reciprocal_h3s_[0] * 60.0},
         {sIndex(1), reciprocal_h3s_[0] * 60.0},
         {vIndex(0), -reciprocal_h2s_[0] * 36.0},
         {vIndex(1), -reciprocal_h2s_[0] * 24.0},
         {aIndex(0), -reciprocal_hs_[0] * 9.0},
         {aIndex(1), reciprocal_hs_[0] * 3.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    qp_problem_.updateConstrainTerm(
        num_samples_ * 4,
        {{sIndex(0), reciprocal_h4s_[0] * 360.0},
         {sIndex(1), -reciprocal_h4s_[0] * 360.0},
         {vIndex(0), reciprocal_h3s_[0] * 192.0},
         {vIndex(1), reciprocal_h3s_[0] * 168.0},
         {aIndex(0), reciprocal_h2s_[0] * 36.0},
         {aIndex(1), -reciprocal_h2s_[0] * 24.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    return;
}

void RouteLineAccModel::setFinalState(double sf, double vf, double af) noexcept {
    if (!std::isnan(sf)) {
        qp_problem_.updateConstrainTerm(
            sIndex(num_samples_ - 1), {{sIndex(num_samples_ - 1), 1.0}},
            sf - tiny_pnc::math::kEpsilon, sf + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            sIndex(num_samples_ - 1), {{sIndex(num_samples_ - 1), 1.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(vf)) {
        qp_problem_.updateConstrainTerm(
            vIndex(num_samples_ - 1), {{vIndex(num_samples_ - 1), 1.0}},
            vf - tiny_pnc::math::kEpsilon, vf + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            vIndex(num_samples_ - 1), {{vIndex(num_samples_ - 1), 1.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(af)) {
        qp_problem_.updateConstrainTerm(
            aIndex(num_samples_ - 1), {{aIndex(num_samples_ - 1), 1.0}},
            af - tiny_pnc::math::kEpsilon, af + tiny_pnc::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            aIndex(num_samples_ - 1), {{aIndex(num_samples_ - 1), 1.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    qp_problem_.updateConstrainTerm(
        num_samples_ * 4 - 1,
        {{sIndex(num_samples_ - 2), -reciprocal_h3s_[num_samples_ - 2] * 60.0},
         {sIndex(num_samples_ - 1), reciprocal_h3s_[num_samples_ - 2] * 60.0},
         {vIndex(num_samples_ - 2), -reciprocal_h2s_[num_samples_ - 2] * 24.0},
         {vIndex(num_samples_ - 1), -reciprocal_h2s_[num_samples_ - 2] * 36.0},
         {aIndex(num_samples_ - 2), -reciprocal_hs_[num_samples_ - 2] * 3.0},
         {aIndex(num_samples_ - 1), reciprocal_hs_[num_samples_ - 2] * 9.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    qp_problem_.updateConstrainTerm(
        num_samples_ * 5 - 1,
        {{sIndex(num_samples_ - 2), -reciprocal_h4s_[num_samples_ - 2] * 360.0},
         {sIndex(num_samples_ - 1), reciprocal_h4s_[num_samples_ - 2] * 360.0},
         {vIndex(num_samples_ - 2), -reciprocal_h3s_[num_samples_ - 2] * 168.0},
         {vIndex(num_samples_ - 1), -reciprocal_h3s_[num_samples_ - 2] * 192.0},
         {aIndex(num_samples_ - 2), -reciprocal_h2s_[num_samples_ - 2] * 24.0},
         {aIndex(num_samples_ - 1), reciprocal_h2s_[num_samples_ - 2] * 36.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    return;
}

void RouteLineAccModel::setVelocityCost(double target_velocity, double velocity_weight) noexcept {
    for (std::size_t i = 0; i < num_samples_ - 1; ++i) {
        const double factor = velocity_weight * reciprocal_hs_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 10.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 20.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 3.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 3.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i), factor * h2s_[i] / 42.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * h2s_[i] / 42.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 10.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 3.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 3.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * h2s_[i] / 42.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * h2s_[i] / 42.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 8.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), -factor * h2s_[i] / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i), factor * h3s_[i] / 30.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i + 1), factor * h3s_[i] / 105.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 8.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i), -factor * h3s_[i] / 105.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * h3s_[i] / 30.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i), factor * h4s_[i] / 630.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i + 1), factor * h4s_[i] / 630.0);
        qp_problem_.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * h4s_[i] / 630.0);
        qp_problem_.addLinCostTerm(sIndex(i), factor * target_velocity * hs_[i] * 2.0);
        qp_problem_.addLinCostTerm(sIndex(i + 1), -factor * target_velocity * hs_[i] * 2.0);
    }
    return;
}

void RouteLineAccModel::setAccelCost(double accel_weight) noexcept {
    for (std::size_t i = 0; i < num_samples_ - 1; ++i) {
        const double factor = accel_weight * reciprocal_h3s_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 240.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i), factor * h2s_[i] * 6.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * h2s_[i] * 6.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 120.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * h2s_[i] * 6.0 / 7.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * h2s_[i] * 6.0 / 7.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 192.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * h2s_[i] * 216.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i), factor * h3s_[i] * 22.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * h3s_[i] * 8.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 192.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * h3s_[i] * 8.0 / 35.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * h3s_[i] * 22.0 / 35.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i), factor * h4s_[i] * 3.0 / 35.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i + 1), factor * h4s_[i] / 35.0);
        qp_problem_.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * h4s_[i] * 3.0 / 35.0);
    }
    return;
}

void RouteLineAccModel::setJerkCost(double jerk_weight) noexcept {
    for (std::size_t i = 0; i < num_samples_ - 1; ++i) {
        const double factor = jerk_weight * reciprocal_h3s_[i] * reciprocal_h2s_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 1440.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i), factor * h2s_[i] * 120.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * h2s_[i] * 120.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 720.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * h2s_[i] * 120.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * h2s_[i] * 120.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 192.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * h2s_[i] * 336.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i), factor * h3s_[i] * 72.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * h3s_[i] * 48.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 192.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * h3s_[i] * 48.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * h3s_[i] * 72.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i), factor * h4s_[i] * 9.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i + 1), -factor * h4s_[i] * 6.0);
        qp_problem_.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * h4s_[i] * 9.0);
    }
    return;
}

void RouteLineAccModel::setSnapCost(double snap_weight) noexcept {
    for (std::size_t i = 0; i < num_samples_ - 1; ++i) {
        const double factor = snap_weight * reciprocal_h4s_[i] * reciprocal_h3s_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 86400.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i), factor * h2s_[i] * 7200.0);
        qp_problem_.addQuadCostTerm(sIndex(i), aIndex(i + 1), -factor * h2s_[i] * 7200.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 43200.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i), -factor * h2s_[i] * 7200.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), aIndex(i + 1), factor * h2s_[i] * 7200.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 10944.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * h2s_[i] * 21312.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i), factor * h3s_[i] * 3744.0);
        qp_problem_.addQuadCostTerm(vIndex(i), aIndex(i + 1), -factor * h3s_[i] * 3456.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 10944.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i), factor * h3s_[i] * 3456.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), aIndex(i + 1), -factor * h3s_[i] * 3744.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i), factor * h4s_[i] * 336.0);
        qp_problem_.addQuadCostTerm(aIndex(i), aIndex(i + 1), -factor * h4s_[i] * 528.0);
        qp_problem_.addQuadCostTerm(aIndex(i + 1), aIndex(i + 1), factor * h4s_[i] * 336.0);
    }
    return;
}

RouteLineAccModel::Result RouteLineAccModel::solve() const noexcept {
    tiny_pnc::math::OsqpSolver osqp_solver{settings_};
    tiny_pnc::math::OsqpResult osqp_result = osqp_solver.solve(qp_problem_);
    const double snap0 =
        (osqp_result.prim_vars[sIndex(0)] - osqp_result.prim_vars[sIndex(1)]) * reciprocal_h4s_[0] *
            360.0 +
        (osqp_result.prim_vars[vIndex(0)] * 192.0 + osqp_result.prim_vars[vIndex(1)] * 168.0) *
            reciprocal_h3s_[0] +
        (osqp_result.prim_vars[aIndex(0)] * 36.0 - osqp_result.prim_vars[aIndex(1)] * 24.0) *
            reciprocal_h2s_[0];
    const double snapf = (osqp_result.prim_vars[sIndex(num_samples_ - 1)] -
                          osqp_result.prim_vars[sIndex(num_samples_ - 2)]) *
                             reciprocal_h4s_[num_samples_ - 2] * 360.0 -
                         (osqp_result.prim_vars[vIndex(num_samples_ - 1)] * 192.0 +
                          osqp_result.prim_vars[vIndex(num_samples_ - 2)] * 168.0) *
                             reciprocal_h3s_[num_samples_ - 2] +
                         (osqp_result.prim_vars[aIndex(num_samples_ - 1)] * 36.0 -
                          osqp_result.prim_vars[aIndex(num_samples_ - 2)] * 24.0) *
                             reciprocal_h2s_[num_samples_ - 2];
    std::array<Motion1d::BoundaryMode, 2> b0{
        Motion1d::BoundaryMode{2, osqp_result.prim_vars[aIndex(0)]},
        Motion1d::BoundaryMode{4, snap0}};
    std::array<Motion1d::BoundaryMode, 2> bf{
        Motion1d::BoundaryMode{2, osqp_result.prim_vars[aIndex(num_samples_ - 1)]},
        Motion1d::BoundaryMode{4, snapf}};
    return Result{
        .motion =
            Motion1d{
                sample_ts_,
                {osqp_result.prim_vars.cbegin(), osqp_result.prim_vars.cbegin() + num_samples_},
                b0,
                bf},
        .info = osqp_result.info};
}

std::size_t RouteLineAccModel::sIndex(std::size_t t_index) const noexcept { return t_index; }

std::size_t RouteLineAccModel::vIndex(std::size_t t_index) const noexcept {
    return num_samples_ + t_index;
}

std::size_t RouteLineAccModel::aIndex(std::size_t t_index) const noexcept {
    return num_samples_ * 2 + t_index;
}

} // namespace kinetics
} // namespace tiny_pnc
