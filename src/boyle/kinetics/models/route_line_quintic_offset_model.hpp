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

#include <array>
#include <format>
#include <limits>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boyle/common/utils/logging.hpp"
#include "boyle/cvxopm/info.hpp"
#include "boyle/cvxopm/problems/qp_problem.hpp"
#include "boyle/cvxopm/settings.hpp"
#include "boyle/cvxopm/solvers/osqp_solver.hpp"
#include "boyle/kinetics/border.hpp"
#include "boyle/kinetics/path.hpp"
#include "boyle/math/curves/piecewise_quintic_curve.hpp"
#include "boyle/math/dense/vec2.hpp"
#include "boyle/math/functions/piecewise_linear_function.hpp"
#include "boyle/math/quintic_interpolation.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::kinetics {

template <::boyle::math::Vec2Arithmetic T>
class RouteLineQuinticOffsetModel final {
  public:
    using value_type = T;
    using scalar_type = typename value_type::value_type;
    using index_type = int;
    using size_type = std::size_t;
    using allocator_type = std::pmr::polymorphic_allocator<value_type>;

    RouteLineQuinticOffsetModel(const RouteLineQuinticOffsetModel& other) noexcept = delete;
    auto operator=(const RouteLineQuinticOffsetModel& other) noexcept
        -> RouteLineQuinticOffsetModel& = delete;
    RouteLineQuinticOffsetModel(RouteLineQuinticOffsetModel&& other) noexcept = delete;
    auto operator=(RouteLineQuinticOffsetModel&& other) noexcept
        -> RouteLineQuinticOffsetModel& = delete;
    ~RouteLineQuinticOffsetModel() noexcept = default;

    [[using gnu: pure, always_inline]]
    auto get_allocator() const noexcept -> allocator_type {
        return m_settings.memory_resource;
    }

    template <
        std::ranges::input_range R0 = std::initializer_list<value_type>,
        std::ranges::input_range R1 = std::initializer_list<scalar_type>>
    explicit RouteLineQuinticOffsetModel(
        R0&& sketch_points, R1&& sample_ss,
        const ::boyle::cvxopm::Settings<scalar_type, index_type>& settings = {}
    ) noexcept(!BOYLE_CHECK_PARAMS)
        requires std::same_as<std::ranges::range_value_t<R0>, value_type> &&
                     std::same_as<std::ranges::range_value_t<R1>, scalar_type>
        : m_num_samples(sample_ss.size()), m_length_scale(sample_ss.back() - sample_ss.front()),
          m_qp_problem(sample_ss.size() * 6, sample_ss.size() * 10, settings.memory_resource),
          m_settings{settings},
          m_sketch_curve(sketch_points, sample_ss.front(), settings.memory_resource),
          m_sample_points(sample_ss.size(), settings.memory_resource),
          m_sample_ss(sample_ss.cbegin(), sample_ss.cend(), settings.memory_resource),
          m_hs(sample_ss.size() - 1, settings.memory_resource),
          m_h2s(sample_ss.size() - 1, settings.memory_resource),
          m_h3s(sample_ss.size() - 1, settings.memory_resource),
          m_h4s(sample_ss.size() - 1, settings.memory_resource),
          m_h6s(sample_ss.size() - 1, settings.memory_resource),
          m_h8s(sample_ss.size() - 1, settings.memory_resource),
          m_reciprocal_hs(sample_ss.size() - 1, settings.memory_resource) {
#if BOYLE_CHECK_PARAMS == 1
        if (sketch_points.size() < 2) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument error detected: size of raw_sketch_points must larger than "
                    "2: "
                    "sketch_points.size() = {0:d}.",
                    sketch_points.size()
                )
            );
        }
        if (sample_ss.size() < 2) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument error detected: size of sample_ss must larger than 2: "
                    "sample_ss.size() = {0:d}.",
                    sample_ss.size()
                )
            );
        }
#endif
        for (size_type i{0}; i < m_sample_ss.size() - 1; ++i) {
            const scalar_type diff = m_sample_ss[i + 1] - m_sample_ss[i];
            const scalar_type reciprocal_diff = 1.0 / diff;
            m_sample_points[i] = m_sketch_curve(m_sample_ss[i]);
            m_hs[i] = diff;
            m_h2s[i] = diff * diff;
            m_h3s[i] = m_h2s[i] * diff;
            m_h4s[i] = m_h2s[i] * m_h2s[i];
            m_h6s[i] = m_h3s[i] * m_h3s[i];
            m_h8s[i] = m_h4s[i] * m_h4s[i];
            m_reciprocal_hs[i] = reciprocal_diff;
        }
        m_sample_points[m_sample_ss.size() - 1] = m_sketch_curve(m_sample_ss.back());
        setIntegrationRelation();
    }

    [[using gnu: pure, always_inline]]
    auto num_samples() const noexcept -> size_type {
        return m_num_samples;
    }

    [[using gnu: pure, always_inline]]
    auto qp_problem() const noexcept -> const ::boyle::cvxopm::QpProblem<scalar_type, index_type>& {
        return m_qp_problem;
    }

    [[using gnu: pure, always_inline]]
    auto settings() const noexcept -> const ::boyle::cvxopm::Settings<scalar_type, index_type>& {
        return m_settings;
    }

    auto setHardBorders(std::span<const HardBorder<value_type>> hard_borders) noexcept -> void {
        std::pmr::vector<scalar_type> x_lower_bound(
            m_num_samples, std::numeric_limits<scalar_type>::lowest(), get_allocator()
        );
        std::pmr::vector<scalar_type> x_upper_bound(
            m_num_samples, std::numeric_limits<scalar_type>::max(), get_allocator()
        );
        std::pmr::vector<scalar_type> y_lower_bound(
            m_num_samples, std::numeric_limits<scalar_type>::lowest(), get_allocator()
        );
        std::pmr::vector<scalar_type> y_upper_bound(
            m_num_samples, std::numeric_limits<scalar_type>::max(), get_allocator()
        );
        for (const HardBorder<value_type>& hard_border : hard_borders) {
            const index_type istart = ::boyle::math::nearestUpperElement(
                                          m_sample_points, hard_border.bound_points.front()
                                      ) -
                                      m_sample_points.cbegin();
            const index_type iend = ::boyle::math::nearestUpperElement(
                                        m_sample_points, hard_border.bound_points.back()
                                    ) -
                                    m_sample_points.cbegin();
            if (istart == m_num_samples) {
                BOYLE_LOG_WARN(
                    "Invalid argument issue detected! This soft border is not in the same region "
                    "as the sketch points."
                );
                continue;
            }
            if (iend == 0) {
                BOYLE_LOG_WARN(
                    "Invalid argument issue detected! This soft border is not in the same region "
                    "as the sketch points."
                );
                continue;
            }

            const auto& bound_points{hard_border.bound_points};
            std::pmr::vector<scalar_type> bound_ss(get_allocator());
            bound_ss.reserve(bound_points.size());
            for (const value_type& bound_point : bound_points) {
                bound_ss.push_back(m_sketch_curve.inverse(bound_point).s);
            }

            const ::boyle::math::pmr::PiecewiseLinearFunction<value_type> bound_func{
                bound_ss, bound_points, get_allocator()
            };

            if (hard_border.chirality == ::boyle::kinetics::Chirality::RIGHT) {
                if (istart != 0) {
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    const scalar_type ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addConstrainTerm(
                        {{xIndex(istart - 1), proration_coeffs[0]},
                         {xIndex(istart), proration_coeffs[1]},
                         {ddxIndex(istart - 1), proration_coeffs[2]},
                         {ddxIndex(istart), proration_coeffs[3]},
                         {d4xIndex(istart - 1), proration_coeffs[4]},
                         {d4xIndex(istart), proration_coeffs[5]}},
                        bound_points.front().x, std::numeric_limits<scalar_type>::max()
                    );
                    m_qp_problem.addConstrainTerm(
                        {{yIndex(istart - 1), proration_coeffs[0]},
                         {yIndex(istart), proration_coeffs[1]},
                         {ddyIndex(istart - 1), proration_coeffs[2]},
                         {ddyIndex(istart), proration_coeffs[3]},
                         {d4xIndex(istart - 1), proration_coeffs[4]},
                         {d4xIndex(istart), proration_coeffs[5]}},
                        bound_points.front().y, std::numeric_limits<scalar_type>::max()
                    );
                }
                for (index_type i{istart}; i < iend; ++i) {
                    const value_type bound_point = bound_func(m_sample_ss[i]);
                    x_lower_bound[i] = std::max(bound_point.x, x_lower_bound[i]);
                    y_lower_bound[i] = std::max(bound_point.y, y_lower_bound[i]);
                }
                if (iend != m_num_samples) {
                    const scalar_type h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                    const scalar_type ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addConstrainTerm(
                        {{xIndex(iend - 1), proration_coeffs[0]},
                         {xIndex(iend), proration_coeffs[1]},
                         {ddxIndex(iend - 1), proration_coeffs[2]},
                         {ddxIndex(iend), proration_coeffs[3]},
                         {d4xIndex(iend - 1), proration_coeffs[4]},
                         {d4xIndex(iend), proration_coeffs[5]}},
                        bound_points.back().x, std::numeric_limits<scalar_type>::max()
                    );
                    m_qp_problem.addConstrainTerm(
                        {{yIndex(iend - 1), proration_coeffs[0]},
                         {yIndex(iend), proration_coeffs[1]},
                         {ddyIndex(iend - 1), proration_coeffs[2]},
                         {ddyIndex(iend), proration_coeffs[3]},
                         {d4yIndex(iend - 1), proration_coeffs[4]},
                         {d4yIndex(iend), proration_coeffs[5]}},
                        bound_points.back().y, std::numeric_limits<scalar_type>::max()
                    );
                }
            } else {
                if (istart != 0) {
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    const scalar_type ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addConstrainTerm(
                        {{xIndex(istart - 1), proration_coeffs[0]},
                         {xIndex(istart), proration_coeffs[1]},
                         {ddxIndex(istart - 1), proration_coeffs[2]},
                         {ddxIndex(istart), proration_coeffs[3]},
                         {d4xIndex(istart - 1), proration_coeffs[4]},
                         {d4xIndex(istart), proration_coeffs[5]}},
                        std::numeric_limits<scalar_type>::lowest(), bound_points.front().x
                    );
                    m_qp_problem.addConstrainTerm(
                        {{yIndex(istart - 1), proration_coeffs[0]},
                         {yIndex(istart), proration_coeffs[1]},
                         {ddyIndex(istart - 1), proration_coeffs[2]},
                         {ddyIndex(istart), proration_coeffs[3]},
                         {d4yIndex(istart - 1), proration_coeffs[4]},
                         {d4yIndex(istart), proration_coeffs[5]}},
                        std::numeric_limits<scalar_type>::lowest(), bound_points.front().y
                    );
                }
                for (index_type i{istart}; i < iend; ++i) {
                    const value_type bound_point = bound_func(m_sample_ss[i]);
                    x_upper_bound[i] = std::min(bound_point.x, x_upper_bound[i]);
                    y_upper_bound[i] = std::min(bound_point.y, y_upper_bound[i]);
                }
                if (iend != m_num_samples) {
                    const scalar_type h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                    const scalar_type ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addConstrainTerm(
                        {{xIndex(iend - 1), proration_coeffs[0]},
                         {xIndex(iend), proration_coeffs[1]},
                         {ddxIndex(iend - 1), proration_coeffs[2]},
                         {ddxIndex(iend), proration_coeffs[3]},
                         {d4xIndex(iend - 1), proration_coeffs[4]},
                         {d4xIndex(iend), proration_coeffs[5]}},
                        std::numeric_limits<scalar_type>::lowest(), bound_points.back().x
                    );
                    m_qp_problem.addConstrainTerm(
                        {{yIndex(iend - 1), proration_coeffs[0]},
                         {yIndex(iend), proration_coeffs[1]},
                         {ddyIndex(iend - 1), proration_coeffs[2]},
                         {ddyIndex(iend), proration_coeffs[3]},
                         {d4yIndex(iend - 1), proration_coeffs[4]},
                         {d4yIndex(iend), proration_coeffs[5]}},
                        std::numeric_limits<scalar_type>::lowest(), bound_points.back().y
                    );
                }
            }
        }

        for (index_type i{1}; i < m_num_samples - 1; ++i) {
            m_qp_problem.updateConstrainTerm(
                xIndex(i), {{xIndex(i), 1.0}}, x_lower_bound[i], x_upper_bound[i]
            );
            m_qp_problem.updateConstrainTerm(
                yIndex(i), {{yIndex(i), 1.0}}, y_lower_bound[i], y_upper_bound[i]
            );
        }
        return;
    }

    auto setSoftBorders(std::span<const SoftBorder<value_type>> soft_borders) noexcept -> void {
        for (const auto& soft_border : soft_borders) {
            const index_type istart = ::boyle::math::nearestUpperElement(
                                          m_sample_points, soft_border.bound_points.front()
                                      ) -
                                      m_sample_points.cbegin();
            const index_type iend = ::boyle::math::nearestUpperElement(
                                        m_sample_points, soft_border.bound_points.back()
                                    ) -
                                    m_sample_points.cbegin();
            if (istart == m_num_samples) {
                BOYLE_LOG_WARN(
                    "Invalid argument issue detected! This soft border is not in the same region "
                    "as the sketch points."
                );
                continue;
            }
            if (iend == 0) {
                BOYLE_LOG_WARN(
                    "Invalid argument issue detected! This soft border is not in the same region "
                    "as the sketch points."
                );
                continue;
            }

            const auto& bound_points{soft_border.bound_points};
            std::pmr::vector<scalar_type> bound_ss(get_allocator());
            bound_ss.reserve(bound_points.size());
            for (const auto& bound_point : bound_points) {
                bound_ss.push_back(m_sketch_curve.inverse(bound_point).s);
            }

            const ::boyle::math::pmr::PiecewiseLinearFunction<value_type> bound_func{
                bound_ss, bound_points, get_allocator()
            };

            scalar_type factor{0.0};

            if (soft_border.chirality == ::boyle::kinetics::Chirality::RIGHT) {
                if (istart == iend) {
                    factor = (bound_ss.back() - bound_ss.front()) / m_length_scale * 0.5;
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    scalar_type ratio{(m_sample_ss.front() - m_sample_ss[istart - 1]) / h};
                    std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(istart - 1), -proration_coeffs[0]},
                         {xIndex(istart), -proration_coeffs[1]},
                         {ddxIndex(istart - 1), -proration_coeffs[2]},
                         {ddxIndex(istart), -proration_coeffs[3]},
                         {d4xIndex(istart - 1), -proration_coeffs[4]},
                         {d4xIndex(istart), -proration_coeffs[5]}},
                        -bound_points.front().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(istart - 1), -proration_coeffs[0]},
                         {yIndex(istart), -proration_coeffs[1]},
                         {ddyIndex(istart - 1), -proration_coeffs[2]},
                         {ddyIndex(istart), -proration_coeffs[3]},
                         {d4yIndex(istart - 1), -proration_coeffs[4]},
                         {d4yIndex(istart), -proration_coeffs[5]}},
                        -bound_points.front().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    ratio = (bound_ss.back() - m_sample_ss[iend - 1]) / h;
                    proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(iend - 1), -proration_coeffs[0]},
                         {xIndex(iend), -proration_coeffs[1]},
                         {ddxIndex(iend - 1), -proration_coeffs[2]},
                         {ddxIndex(iend), -proration_coeffs[3]},
                         {d4xIndex(iend - 1), -proration_coeffs[4]},
                         {d4xIndex(iend), -proration_coeffs[5]}},
                        -bound_points.back().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(iend - 1), -proration_coeffs[0]},
                         {yIndex(iend), -proration_coeffs[1]},
                         {ddyIndex(iend - 1), -proration_coeffs[2]},
                         {ddyIndex(iend), -proration_coeffs[3]},
                         {d4yIndex(iend - 1), -proration_coeffs[4]},
                         {d4yIndex(iend), -proration_coeffs[5]}},
                        -bound_points.back().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    continue;
                }
                if (istart != 0) {
                    factor = (m_sample_ss[istart] - bound_ss.front()) / m_length_scale;
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    const scalar_type ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(istart - 1), -proration_coeffs[0]},
                         {xIndex(istart), -proration_coeffs[1]},
                         {ddxIndex(istart - 1), -proration_coeffs[2]},
                         {ddxIndex(istart), -proration_coeffs[3]},
                         {d4xIndex(istart - 1), -proration_coeffs[4]},
                         {d4xIndex(istart), -proration_coeffs[5]}},
                        -bound_points.front().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(istart - 1), -proration_coeffs[0]},
                         {yIndex(istart), -proration_coeffs[1]},
                         {ddyIndex(istart - 1), -proration_coeffs[2]},
                         {ddyIndex(istart), -proration_coeffs[3]},
                         {d4yIndex(istart - 1), -proration_coeffs[4]},
                         {d4yIndex(istart), -proration_coeffs[5]}},
                        -bound_points.front().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
                for (index_type i{istart}; i < iend; ++i) {
                    if (i == istart) {
                        factor = m_hs[i] / m_length_scale * 0.5;
                    } else if (i == iend - 1) {
                        factor = m_hs[i - 1] / m_length_scale * 0.5;
                    } else {
                        factor = (m_hs[i - 1] + m_hs[i]) / m_length_scale * 0.5;
                    }
                    const value_type bound_point = bound_func(m_sample_ss[i]);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(i), -1.0}}, -bound_point.x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(i), -1.0}}, -bound_point.y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
                if (iend != m_num_samples) {
                    factor = (bound_ss.back() - m_sample_ss[iend - 1]) / m_length_scale;
                    const scalar_type h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                    const scalar_type ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(iend - 1), -proration_coeffs[0]},
                         {xIndex(iend), -proration_coeffs[1]},
                         {ddxIndex(iend - 1), -proration_coeffs[2]},
                         {ddxIndex(iend), -proration_coeffs[3]},
                         {d4xIndex(iend - 1), -proration_coeffs[4]},
                         {d4xIndex(iend), -proration_coeffs[5]}},
                        -bound_points.back().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(iend - 1), -proration_coeffs[0]},
                         {yIndex(iend), -proration_coeffs[1]},
                         {ddyIndex(iend - 1), -proration_coeffs[2]},
                         {ddyIndex(iend), -proration_coeffs[3]},
                         {d4yIndex(iend - 1), -proration_coeffs[4]},
                         {d4yIndex(iend), -proration_coeffs[5]}},
                        -bound_points.back().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
            } else {
                if (istart == iend) {
                    factor = (bound_ss.back() - bound_ss.front()) / m_length_scale * 0.5;
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    scalar_type ratio{(m_sample_ss.front() - m_sample_ss[istart - 1]) / h};
                    std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(istart - 1), proration_coeffs[0]},
                         {xIndex(istart), proration_coeffs[1]},
                         {ddxIndex(istart - 1), proration_coeffs[2]},
                         {ddxIndex(istart), proration_coeffs[3]},
                         {d4xIndex(istart - 1), proration_coeffs[4]},
                         {d4xIndex(istart), proration_coeffs[5]}},
                        bound_points.front().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(istart - 1), proration_coeffs[0]},
                         {yIndex(istart), proration_coeffs[1]},
                         {ddyIndex(istart - 1), proration_coeffs[2]},
                         {ddyIndex(istart), proration_coeffs[3]},
                         {d4yIndex(istart - 1), proration_coeffs[4]},
                         {d4yIndex(istart), proration_coeffs[5]}},
                        bound_points.front().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    ratio = (bound_ss.back() - m_sample_ss[iend - 1]) / h;
                    proration_coeffs = ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(iend - 1), proration_coeffs[0]},
                         {xIndex(iend), proration_coeffs[1]},
                         {ddxIndex(iend - 1), proration_coeffs[2]},
                         {ddxIndex(iend), proration_coeffs[3]},
                         {d4xIndex(iend - 1), proration_coeffs[4]},
                         {d4xIndex(iend), proration_coeffs[5]}},
                        bound_points.back().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(iend - 1), proration_coeffs[0]},
                         {yIndex(iend), proration_coeffs[1]},
                         {ddyIndex(iend - 1), proration_coeffs[2]},
                         {ddyIndex(iend), proration_coeffs[3]},
                         {d4yIndex(iend - 1), proration_coeffs[4]},
                         {d4yIndex(iend), proration_coeffs[5]}},
                        bound_points.back().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    continue;
                }
                if (istart != 0) {
                    factor = (m_sample_ss[istart] - bound_ss.front()) / m_length_scale;
                    const scalar_type h{m_sample_ss[istart] - m_sample_ss[istart - 1]};
                    const scalar_type ratio{(bound_ss.front() - m_sample_ss[istart - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(istart - 1), proration_coeffs[0]},
                         {xIndex(istart), proration_coeffs[1]},
                         {ddxIndex(istart - 1), proration_coeffs[2]},
                         {ddxIndex(istart), proration_coeffs[3]},
                         {d4xIndex(istart - 1), proration_coeffs[4]},
                         {d4xIndex(istart), proration_coeffs[5]}},
                        bound_points.front().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(istart - 1), proration_coeffs[0]},
                         {yIndex(istart), proration_coeffs[1]},
                         {ddyIndex(istart - 1), proration_coeffs[2]},
                         {ddyIndex(istart), proration_coeffs[3]},
                         {d4yIndex(istart - 1), proration_coeffs[4]},
                         {d4yIndex(istart), proration_coeffs[5]}},
                        bound_points.front().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
                for (index_type i{istart}; i < iend; ++i) {
                    if (i == istart) {
                        factor = m_hs[i] / m_length_scale * 0.5;
                    } else if (i == iend - 1) {
                        factor = m_hs[i - 1] / m_length_scale * 0.5;
                    } else {
                        factor = (m_hs[i - 1] + m_hs[i]) / m_length_scale * 0.5;
                    }
                    const value_type bound_point = bound_func(m_sample_ss[i]);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(i), 1.0}}, bound_point.x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(i), 1.0}}, bound_point.y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
                if (iend != m_num_samples) {
                    factor = (bound_ss.back() - m_sample_ss[iend - 1]) / m_length_scale;
                    const scalar_type h{m_sample_ss[iend] - m_sample_ss[iend - 1]};
                    const scalar_type ratio{(bound_ss.back() - m_sample_ss[iend - 1]) / h};
                    const std::array<scalar_type, 6> proration_coeffs =
                        ::boyle::math::quinerpCoeffs(ratio, h);
                    m_qp_problem.addRampCostTerm(
                        {{xIndex(iend - 1), proration_coeffs[0]},
                         {xIndex(iend), proration_coeffs[1]},
                         {ddxIndex(iend - 1), proration_coeffs[2]},
                         {ddxIndex(iend), proration_coeffs[3]},
                         {d4xIndex(iend - 1), proration_coeffs[4]},
                         {d4xIndex(iend), proration_coeffs[5]}},
                        bound_points.back().x, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                    m_qp_problem.addRampCostTerm(
                        {{yIndex(iend - 1), proration_coeffs[0]},
                         {yIndex(iend), proration_coeffs[1]},
                         {ddyIndex(iend - 1), proration_coeffs[2]},
                         {ddyIndex(iend), proration_coeffs[3]},
                         {d4yIndex(iend - 1), proration_coeffs[4]},
                         {d4yIndex(iend), proration_coeffs[5]}},
                        bound_points.back().y, soft_border.linear_weight * factor,
                        soft_border.quadratic_weight * factor
                    );
                }
            }
        }
        return;
    }

    auto setDdxRange(scalar_type ddx_min, scalar_type ddx_max) noexcept -> void {
        for (index_type i{0}; i < m_num_samples; i++) {
            m_qp_problem.updateConstrainTerm(ddxIndex(i), {{ddxIndex(i), 1.0}}, ddx_min, ddx_max);
        }
        return;
    }

    auto setDdyRange(scalar_type ddy_min, scalar_type ddy_max) noexcept -> void {
        for (index_type i{0}; i < m_num_samples; i++) {
            m_qp_problem.updateConstrainTerm(ddyIndex(i), {{ddyIndex(i), 1.0}}, ddy_min, ddy_max);
        }
        return;
    }

    auto setD4xRange(scalar_type d4x_min, scalar_type d4x_max) noexcept -> void {
        for (index_type i{0}; i < m_num_samples; i++) {
            m_qp_problem.updateConstrainTerm(d4xIndex(i), {{d4xIndex(i), 1.0}}, d4x_min, d4x_max);
        }
        return;
    }

    auto setD4yRange(scalar_type d4y_min, scalar_type d4y_max) noexcept -> void {
        for (index_type i{0}; i < m_num_samples; i++) {
            m_qp_problem.updateConstrainTerm(d4yIndex(i), {{d4yIndex(i), 1.0}}, d4y_min, d4y_max);
        }
        return;
    }

    auto setInitialState(
        value_type r0,
        value_type t0 =
            {std::numeric_limits<scalar_type>::quiet_NaN(),
             std::numeric_limits<scalar_type>::quiet_NaN()},
        value_type n0 =
            {std::numeric_limits<scalar_type>::quiet_NaN(),
             std::numeric_limits<scalar_type>::quiet_NaN()},
        value_type j0 = {
            std::numeric_limits<scalar_type>::quiet_NaN(),
            std::numeric_limits<scalar_type>::quiet_NaN()
        }
    ) noexcept -> void {
        m_qp_problem.updateConstrainTerm(xIndex(0), {{xIndex(0), 1.0}}, r0.x, r0.x);
        if (!std::isnan(t0.x)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 6,
                {{xIndex(0), -m_reciprocal_hs[0]},
                 {xIndex(1), m_reciprocal_hs[0]},
                 {ddxIndex(0), -m_hs[0] / 3.0},
                 {ddxIndex(1), -m_hs[0] / 6.0},
                 {d4xIndex(0), m_h3s[0] / 45.0},
                 {d4xIndex(1), m_h3s[0] * 7.0 / 360.0}},
                t0.x, t0.x
            );
        }
        if (!std::isnan(n0.x)) {
            m_qp_problem.updateConstrainTerm(ddxIndex(0), {{ddxIndex(0), 1.0}}, n0.x, n0.x);
        }
        if (!std::isnan(j0.x)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 7,
                {{ddxIndex(0), -m_reciprocal_hs[0]},
                 {ddxIndex(1), m_reciprocal_hs[0]},
                 {d4xIndex(0), -m_hs[0] / 3.0},
                 {d4xIndex(1), -m_hs[0] / 6.0}},
                j0.x, j0.x
            );
        }
        // m_qp_problem.updateConstrainTerm(
        //     d4xIndex(0), {{d4xIndex(0), 1.0}}, std::numeric_limits<scalar_type>::lowest(),
        //     std::numeric_limits<scalar_type>::max()
        // );

        m_qp_problem.updateConstrainTerm(yIndex(0), {{yIndex(0), 1.0}}, r0.y, r0.y);
        if (!std::isnan(t0.y)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 8,
                {{yIndex(0), -m_reciprocal_hs[0]},
                 {yIndex(1), m_reciprocal_hs[0]},
                 {ddyIndex(0), -m_hs[0] / 3.0},
                 {ddyIndex(1), -m_hs[0] / 6.0},
                 {d4yIndex(0), m_h3s[0] / 45.0},
                 {d4yIndex(1), m_h3s[0] * 7.0 / 360.0}},
                t0.y, t0.y
            );
        }
        if (!std::isnan(n0.y)) {
            m_qp_problem.updateConstrainTerm(ddyIndex(0), {{ddyIndex(0), 1.0}}, n0.y, n0.y);
        }
        if (!std::isnan(j0.y)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 9,
                {{ddyIndex(0), -m_reciprocal_hs[0]},
                 {ddyIndex(1), m_reciprocal_hs[0]},
                 {d4yIndex(0), -m_hs[0] / 3.0},
                 {d4yIndex(1), -m_hs[0] / 6.0}},
                j0.y, j0.y
            );
        }
        // m_qp_problem.updateConstrainTerm(
        //     d4yIndex(0), {{d4yIndex(0), 1.0}}, std::numeric_limits<scalar_type>::lowest(),
        //     std::numeric_limits<scalar_type>::max()
        // );
        return;
    }

    auto setFinalState(
        value_type rf,
        value_type tf =
            {std::numeric_limits<scalar_type>::quiet_NaN(),
             std::numeric_limits<scalar_type>::quiet_NaN()},
        value_type nf =
            {std::numeric_limits<scalar_type>::quiet_NaN(),
             std::numeric_limits<scalar_type>::quiet_NaN()},
        value_type jf = {
            std::numeric_limits<scalar_type>::quiet_NaN(),
            std::numeric_limits<scalar_type>::quiet_NaN()
        }
    ) noexcept -> void {
        m_qp_problem.updateConstrainTerm(
            xIndex(m_num_samples - 1), {{xIndex(m_num_samples - 1), 1.0}}, rf.x, rf.x
        );
        if (!std::isnan(tf.x)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 7 - 1,
                {{xIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
                 {xIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
                 {ddxIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
                 {ddxIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0},
                 {d4xIndex(m_num_samples - 2), -m_h3s[m_num_samples - 2] * 7.0 / 360.0},
                 {d4xIndex(m_num_samples - 1), -m_h3s[m_num_samples - 2] / 45.0}},
                tf.x, tf.x
            );
        }
        if (!std::isnan(nf.x)) {
            m_qp_problem.updateConstrainTerm(
                ddxIndex(m_num_samples - 1), {{ddxIndex(m_num_samples - 1), 1.0}}, nf.x, nf.x
            );
        }
        if (!std::isnan(jf.x)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 8 - 1,
                {{ddxIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
                 {ddxIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
                 {d4xIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
                 {d4xIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0}},
                jf.x, jf.x
            );
        }
        // m_qp_problem.updateConstrainTerm(
        //     d4xIndex(m_num_samples - 1), {{d4xIndex(m_num_samples - 1), 1.0}},
        //     std::numeric_limits<scalar_type>::lowest(), std::numeric_limits<scalar_type>::max()
        // );

        m_qp_problem.updateConstrainTerm(
            yIndex(m_num_samples - 1), {{yIndex(m_num_samples - 1), 1.0}}, rf.y, rf.y
        );
        if (!std::isnan(tf.y)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 9 - 1,
                {{yIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
                 {yIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
                 {ddyIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
                 {ddyIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0},
                 {d4yIndex(m_num_samples - 2), -m_h3s[m_num_samples - 2] * 7.0 / 360.0},
                 {d4yIndex(m_num_samples - 1), -m_h3s[m_num_samples - 2] / 45.0}},
                tf.y, tf.y
            );
        }
        if (!std::isnan(nf.y)) {
            m_qp_problem.updateConstrainTerm(
                ddyIndex(m_num_samples - 1), {{ddyIndex(m_num_samples - 1), 1.0}}, nf.y, nf.y
            );
        }
        if (!std::isnan(jf.y)) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 8 - 1,
                {{ddyIndex(m_num_samples - 2), -m_reciprocal_hs[m_num_samples - 2]},
                 {ddyIndex(m_num_samples - 1), m_reciprocal_hs[m_num_samples - 2]},
                 {d4yIndex(m_num_samples - 2), m_hs[m_num_samples - 2] / 6.0},
                 {d4yIndex(m_num_samples - 1), m_hs[m_num_samples - 2] / 3.0}},
                jf.y, jf.y
            );
        }
        // m_qp_problem.updateConstrainTerm(
        //     d4yIndex(m_num_samples - 1), {{d4yIndex(m_num_samples - 1), 1.0}},
        //     std::numeric_limits<scalar_type>::lowest(), std::numeric_limits<scalar_type>::max()
        // );
        return;
    }

    auto setOffsetCost(scalar_type offset_weight) noexcept -> void {
        constexpr std::array<scalar_type, 10> kFactors{
            1.0 / 3.0,     2.0 / 45.0,   7.0 / 180.0,      2.0 / 945.0,   4.0 / 945.0,
            31.0 / 7560.0, 2.0 / 4725.0, 127.0 / 302400.0, 2.0 / 93555.0, 73.0 / 1710720.0
        };
        for (index_type i{0}; i < m_num_samples - 1; ++i) {
            const scalar_type factor = offset_weight * m_hs[i] / m_length_scale;

            m_qp_problem.addQuadCostTerm(xIndex(i), xIndex(i), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(xIndex(i), xIndex(i + 1), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(xIndex(i), ddxIndex(i), -factor * m_h2s[i] * kFactors[1]);
            m_qp_problem.addQuadCostTerm(
                xIndex(i), ddxIndex(i + 1), -factor * m_h2s[i] * kFactors[2]
            );
            m_qp_problem.addQuadCostTerm(xIndex(i), d4xIndex(i), factor * m_h4s[i] * kFactors[4]);
            m_qp_problem.addQuadCostTerm(
                xIndex(i), d4xIndex(i + 1), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(xIndex(i + 1), xIndex(i + 1), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(
                xIndex(i + 1), ddxIndex(i), -factor * m_h2s[i] * kFactors[2]
            );
            m_qp_problem.addQuadCostTerm(
                xIndex(i + 1), ddxIndex(i + 1), -factor * m_h2s[i] * kFactors[1]
            );
            m_qp_problem.addQuadCostTerm(
                xIndex(i + 1), d4xIndex(i), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(
                xIndex(i + 1), d4xIndex(i + 1), factor * m_h4s[i] * kFactors[4]
            );
            m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i] * kFactors[3]);
            m_qp_problem.addQuadCostTerm(
                ddxIndex(i), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(
                ddxIndex(i), d4xIndex(i), -factor * m_h6s[i] * kFactors[6]
            );
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
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[9]
            );
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[8]
            );
            m_qp_problem.addLinCostTerm(
                xIndex(i),
                -factor * (m_sample_points[i].x * 2.0 + m_sample_points[i + 1].x) * kFactors[0]
            );
            m_qp_problem.addLinCostTerm(
                xIndex(i + 1),
                -factor * (m_sample_points[i].x + m_sample_points[i + 1].x * 2.0) * kFactors[0]
            );
            m_qp_problem.addLinCostTerm(
                ddxIndex(i),
                factor * m_h2s[i] *
                    (m_sample_points[i].x * kFactors[1] + m_sample_points[i + 1].x * kFactors[2])
            );
            m_qp_problem.addLinCostTerm(
                ddxIndex(i + 1),
                factor * m_h2s[i] *
                    (m_sample_points[i].x * kFactors[2] + m_sample_points[i + 1].x * kFactors[1])
            );
            m_qp_problem.addLinCostTerm(
                d4xIndex(i),
                -factor * m_h4s[i] *
                    (m_sample_points[i].x * kFactors[4] + m_sample_points[i + 1].x * kFactors[5])
            );
            m_qp_problem.addLinCostTerm(
                d4xIndex(i + 1),
                -factor * m_h4s[i] *
                    (m_sample_points[i].x * kFactors[5] + m_sample_points[i + 1].x * kFactors[4])
            );

            m_qp_problem.addQuadCostTerm(yIndex(i), yIndex(i), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(yIndex(i), yIndex(i + 1), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(yIndex(i), ddyIndex(i), -factor * m_h2s[i] * kFactors[1]);
            m_qp_problem.addQuadCostTerm(
                yIndex(i), ddyIndex(i + 1), -factor * m_h2s[i] * kFactors[2]
            );
            m_qp_problem.addQuadCostTerm(yIndex(i), d4yIndex(i), factor * m_h4s[i] * kFactors[4]);
            m_qp_problem.addQuadCostTerm(
                yIndex(i), d4yIndex(i + 1), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(yIndex(i + 1), yIndex(i + 1), factor * kFactors[0]);
            m_qp_problem.addQuadCostTerm(
                yIndex(i + 1), ddyIndex(i), -factor * m_h2s[i] * kFactors[2]
            );
            m_qp_problem.addQuadCostTerm(
                yIndex(i + 1), ddyIndex(i + 1), -factor * m_h2s[i] * kFactors[1]
            );
            m_qp_problem.addQuadCostTerm(
                yIndex(i + 1), d4yIndex(i), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(
                yIndex(i + 1), d4yIndex(i + 1), factor * m_h4s[i] * kFactors[4]
            );
            m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i] * kFactors[3]);
            m_qp_problem.addQuadCostTerm(
                ddyIndex(i), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[5]
            );
            m_qp_problem.addQuadCostTerm(
                ddyIndex(i), d4yIndex(i), -factor * m_h6s[i] * kFactors[6]
            );
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
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[9]
            );
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[8]
            );
            m_qp_problem.addLinCostTerm(
                yIndex(i),
                -factor * (m_sample_points[i].y * 2.0 + m_sample_points[i + 1].y) * kFactors[0]
            );
            m_qp_problem.addLinCostTerm(
                yIndex(i + 1),
                -factor * (m_sample_points[i].y + m_sample_points[i + 1].y * 2.0) * kFactors[0]
            );
            m_qp_problem.addLinCostTerm(
                ddyIndex(i),
                factor * m_h2s[i] *
                    (m_sample_points[i].y * kFactors[1] + m_sample_points[i + 1].y * kFactors[2])
            );
            m_qp_problem.addLinCostTerm(
                ddyIndex(i + 1),
                factor * m_h2s[i] *
                    (m_sample_points[i].y * kFactors[2] + m_sample_points[i + 1].y * kFactors[1])
            );
            m_qp_problem.addLinCostTerm(
                d4yIndex(i),
                -factor * m_h4s[i] *
                    (m_sample_points[i].y * kFactors[4] + m_sample_points[i + 1].y * kFactors[5])
            );
            m_qp_problem.addLinCostTerm(
                d4yIndex(i + 1),
                -factor * m_h4s[i] *
                    (m_sample_points[i].y * kFactors[5] + m_sample_points[i + 1].y * kFactors[4])
            );
        }
        return;
    }

    auto setCurvatureCost(scalar_type curvature_weight) noexcept -> void {
        constexpr std::array<scalar_type, 5> kFactors{
            1.0 / 3.0, 2.0 / 45.0, 7.0 / 180.0, 2.0 / 945.0, 31.0 / 7560.0
        };
        for (index_type i{0}; i < m_num_samples - 1; ++i) {
            const scalar_type factor = curvature_weight * m_reciprocal_hs[i] * m_reciprocal_hs[i] *
                                       m_reciprocal_hs[i] / m_length_scale;

            m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i] * kFactors[0]);
            m_qp_problem.addQuadCostTerm(
                ddxIndex(i), ddxIndex(i + 1), factor * m_h4s[i] * kFactors[0]
            );
            m_qp_problem.addQuadCostTerm(
                ddxIndex(i), d4xIndex(i), -factor * m_h6s[i] * kFactors[1]
            );
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
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[4]
            );
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] * kFactors[3]
            );

            m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i] * kFactors[0]);
            m_qp_problem.addQuadCostTerm(
                ddyIndex(i), ddyIndex(i + 1), factor * m_h4s[i] * kFactors[0]
            );
            m_qp_problem.addQuadCostTerm(
                ddyIndex(i), d4yIndex(i), -factor * m_h6s[i] * kFactors[1]
            );
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
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[4]
            );
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] * kFactors[3]
            );
        }
        return;
    }

    auto setDCurvatureCost(scalar_type dcurvature_weight) noexcept -> void {
        for (index_type i{0}; i < m_num_samples - 1; ++i) {
            const scalar_type factor = dcurvature_weight * m_reciprocal_hs[i] * m_reciprocal_hs[i] *
                                       m_reciprocal_hs[i] * m_reciprocal_hs[i] *
                                       m_reciprocal_hs[i] / m_length_scale;
            m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i), factor * m_h4s[i]);
            m_qp_problem.addQuadCostTerm(ddxIndex(i), ddxIndex(i + 1), -factor * m_h4s[i] * 2.0);
            m_qp_problem.addQuadCostTerm(ddxIndex(i + 1), ddxIndex(i + 1), factor * m_h4s[i]);
            m_qp_problem.addQuadCostTerm(d4xIndex(i), d4xIndex(i), factor * m_h8s[i] / 45.0);
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i), d4xIndex(i + 1), factor * m_h8s[i] * 7.0 / 180.0
            );
            m_qp_problem.addQuadCostTerm(
                d4xIndex(i + 1), d4xIndex(i + 1), factor * m_h8s[i] / 45.0
            );

            m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i), factor * m_h4s[i]);
            m_qp_problem.addQuadCostTerm(ddyIndex(i), ddyIndex(i + 1), -factor * m_h4s[i] * 2.0);
            m_qp_problem.addQuadCostTerm(ddyIndex(i + 1), ddyIndex(i + 1), factor * m_h4s[i]);
            m_qp_problem.addQuadCostTerm(d4yIndex(i), d4yIndex(i), factor * m_h8s[i] / 45.0);
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i), d4yIndex(i + 1), factor * m_h8s[i] * 7.0 / 180.0
            );
            m_qp_problem.addQuadCostTerm(
                d4yIndex(i + 1), d4yIndex(i + 1), factor * m_h8s[i] / 45.0
            );
        }
        return;
    }

    auto solve() const
        -> std::pair<Path<value_type>, ::boyle::cvxopm::Info<scalar_type, index_type>> {
        const ::boyle::cvxopm::OsqpSolver<scalar_type, index_type> osqp_solver{m_settings};
        const auto [osqp_result, osqp_info] = osqp_solver.solve(m_qp_problem);
        const scalar_type ddx0{osqp_result.prim_vars[ddxIndex(0)]};
        const scalar_type ddxf{osqp_result.prim_vars[ddxIndex(m_num_samples - 1)]};
        const scalar_type ddy0{osqp_result.prim_vars[ddyIndex(0)]};
        const scalar_type ddyf{osqp_result.prim_vars[ddyIndex(m_num_samples - 1)]};
        const scalar_type d4x0{osqp_result.prim_vars[d4xIndex(0)]};
        const scalar_type d4xf{osqp_result.prim_vars[d4xIndex(m_num_samples - 1)]};
        const scalar_type d4y0{osqp_result.prim_vars[d4yIndex(0)]};
        const scalar_type d4yf{osqp_result.prim_vars[d4yIndex(m_num_samples - 1)]};
        const std::array<::boyle::math::BoundaryMode<value_type>, 2> b0{
            ::boyle::math::BoundaryMode<value_type>{2, {ddx0, ddy0}},
            ::boyle::math::BoundaryMode<value_type>{4, {d4x0, d4y0}}
        };
        const std::array<::boyle::math::BoundaryMode<value_type>, 2> bf{
            ::boyle::math::BoundaryMode<value_type>{2, {ddxf, ddyf}},
            ::boyle::math::BoundaryMode<value_type>{4, {d4xf, d4yf}}
        };
        std::pmr::vector<value_type> anchor_points(m_num_samples, get_allocator());
        ::boyle::math::squeeze(
            osqp_result.prim_vars.data() + xIndex(0),
            osqp_result.prim_vars.data() + xIndex(m_num_samples - 1) + 1,
            osqp_result.prim_vars.data() + yIndex(0),
            osqp_result.prim_vars.data() + yIndex(m_num_samples - 1) + 1, anchor_points.begin()
        );
        return std::make_pair(
            Path<value_type>{anchor_points, b0, bf, m_sample_ss.front(), get_allocator()}, osqp_info
        );
    }

  private:
    auto setIntegrationRelation() noexcept -> void {
        for (index_type i{1}; i < m_num_samples - 1; ++i) {
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 6 + i,
                {{xIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
                 {xIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
                 {xIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
                 {ddxIndex(i - 1), m_hs[i - 1]},
                 {ddxIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
                 {ddxIndex(i + 1), m_hs[i]},
                 {d4xIndex(i - 1), -m_h3s[i - 1] * 7.0 / 60.0},
                 {d4xIndex(i), -(m_h3s[i - 1] + m_h3s[i]) * 2.0 / 15.0},
                 {d4xIndex(i + 1), -m_h3s[i] * 7.0 / 60.0}},
                -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
            );
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 7 + i,
                {{ddxIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
                 {ddxIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
                 {ddxIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
                 {d4xIndex(i - 1), m_hs[i - 1]},
                 {d4xIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
                 {d4xIndex(i + 1), m_hs[i]}},
                -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
            );
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 8 + i,
                {{yIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
                 {yIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
                 {yIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
                 {ddyIndex(i - 1), m_hs[i - 1]},
                 {ddyIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
                 {ddyIndex(i + 1), m_hs[i]},
                 {d4yIndex(i - 1), -m_h3s[i - 1] * 7.0 / 60.0},
                 {d4yIndex(i), -(m_h3s[i - 1] + m_h3s[i]) * 2.0 / 15.0},
                 {d4yIndex(i + 1), -m_h3s[i] * 7.0 / 60.0}},
                -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
            );
            m_qp_problem.updateConstrainTerm(
                m_num_samples * 9 + i,
                {{ddyIndex(i - 1), -m_reciprocal_hs[i - 1] * 6.0},
                 {ddyIndex(i), (m_reciprocal_hs[i - 1] + m_reciprocal_hs[i]) * 6.0},
                 {ddyIndex(i + 1), -m_reciprocal_hs[i] * 6.0},
                 {d4yIndex(i - 1), m_hs[i - 1]},
                 {d4yIndex(i), (m_hs[i - 1] + m_hs[i]) * 2.0},
                 {d4yIndex(i + 1), m_hs[i]}},
                -::boyle::math::kEpsilon, ::boyle::math::kEpsilon
            );
        }
        return;
    }

    [[using gnu: pure, always_inline]]
    auto xIndex(index_type s_index) const noexcept -> index_type {
        return s_index;
    }

    [[using gnu: pure, always_inline]]
    auto ddxIndex(index_type s_index) const noexcept -> index_type {
        return m_num_samples + s_index;
    }

    [[using gnu: pure, always_inline]]
    auto d4xIndex(index_type s_index) const noexcept -> index_type {
        return m_num_samples * 2 + s_index;
    }

    [[using gnu: pure, always_inline]]
    auto yIndex(index_type s_index) const noexcept -> index_type {
        return m_num_samples * 3 + s_index;
    }

    [[using gnu: pure, always_inline]]
    auto ddyIndex(index_type s_index) const noexcept -> index_type {
        return m_num_samples * 4 + s_index;
    }

    [[using gnu: pure, always_inline]]
    auto d4yIndex(index_type s_index) const noexcept -> index_type {
        return m_num_samples * 5 + s_index;
    }

    index_type m_num_samples;
    scalar_type m_length_scale;
    ::boyle::cvxopm::QpProblem<scalar_type, index_type> m_qp_problem;
    ::boyle::cvxopm::Settings<scalar_type, index_type> m_settings;
    ::boyle::math::pmr::PiecewiseQuinticCurve<value_type> m_sketch_curve;
    std::pmr::vector<value_type> m_sample_points;
    std::pmr::vector<scalar_type> m_sample_ss;
    std::pmr::vector<scalar_type> m_hs;
    std::pmr::vector<scalar_type> m_h2s;
    std::pmr::vector<scalar_type> m_h3s;
    std::pmr::vector<scalar_type> m_h4s;
    std::pmr::vector<scalar_type> m_h6s;
    std::pmr::vector<scalar_type> m_h8s;
    std::pmr::vector<scalar_type> m_reciprocal_hs;
};

} // namespace boyle::kinetics
