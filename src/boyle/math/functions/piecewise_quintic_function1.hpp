/**
 * @file piecewise_quintic_function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <format>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/math/concepts.hpp"
#include "boyle/math/dense/vec2.hpp"
#include "boyle/math/dense/vec3.hpp"
#include "boyle/math/quintic_interpolation.hpp"
#include "boyle/math/utils.hpp"

namespace boyle::math {

template <GeneralArithmetic T, std::floating_point U = T>
class PiecewiseQuinticFunction1 final {
    friend class boost::serialization::access;

  public:
    using value_type = T;
    using param_type = U;

    static constexpr param_type kDuplicateCriterion{1E-8};

    struct [[nodiscard]] BoundaryMode final {
        unsigned int order;
        value_type derivative;
    };

    PiecewiseQuinticFunction1() noexcept = default;
    PiecewiseQuinticFunction1(const PiecewiseQuinticFunction1& other) noexcept = default;
    auto operator=(const PiecewiseQuinticFunction1& other) noexcept
        -> PiecewiseQuinticFunction1& = default;
    PiecewiseQuinticFunction1(PiecewiseQuinticFunction1&& other) noexcept = default;
    auto operator=(PiecewiseQuinticFunction1&& other) noexcept
        -> PiecewiseQuinticFunction1& = default;
    ~PiecewiseQuinticFunction1() noexcept = default;

    [[using gnu: always_inline]]
    explicit PiecewiseQuinticFunction1(std::vector<param_type> ts, std::vector<value_type> ys)
        : PiecewiseQuinticFunction1(
              std::move(ts), std::move(ys),
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}},
              std::array<BoundaryMode, 2>{BoundaryMode{2, T{0.0}}, BoundaryMode{4, T{0.0}}}
          ) {}

    [[using gnu: ]]
    explicit PiecewiseQuinticFunction1(
        std::vector<param_type> ts, std::vector<value_type> ys, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf
    ) noexcept(!BOYLE_CHECK_PARAMS)
        : m_ts{std::move(ts)}, m_ys{std::move(ys)} {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ts.size() < 2 || m_ys.size() < 2) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() "
                    "= "
                    "{0:d} while ys.size() = {1:d}",
                    m_ts.size(), m_ys.size()
                )
            );
        }
        if (m_ts.size() != m_ys.size()) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! ts, ys must share the same size: ts.size() = "
                    "{0:d} "
                    "while ys.size() = {1:d}",
                    m_ts.size(), m_ys.size()
                )
            );
        }
        if (!std::is_sorted(m_ts.cbegin(), m_ts.cend())) [[unlikely]] {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        }
        if (hasDuplicates(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, kDuplicateCriterion))
            [[unlikely]] {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }

        if (b0[0].order == 0 || b0[0].order > 4) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! The derivative order of b0 can only be either 1, "
                    "2, 3, "
                    "4, 5: b01.order = {0:d}.",
                    b0[0].order
                )
            );
        }
        if (b0[1].order == 0 || b0[1].order > 4) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! The derivative order of b0 can only be either 1, "
                    "2, 3, "
                    "4, 5: b02.order = {0:d}.",
                    b0[0].order
                )
            );
        }
        if (b0[0].order == b0[1].order) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! You should not set the same order of boundary mode "
                    "twice: b01.order == {0:d} while b02.order == {1:d}",
                    b0[0].order, b0[1].order
                )
            );
        }
        if (bf[0].order == 0 || bf[0].order > 4) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! The derivative order of b0 can only be either 1, "
                    "2, 3, "
                    "4, 5: bf1.order = {0:d}.",
                    bf[0].order
                )
            );
        }
        if (bf[1].order == 0 || bf[1].order > 4) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! The derivative order of b0 can only be either 1, "
                    "2, 3, "
                    "4, 5: bf2.order = {0:d}.",
                    bf[1].order
                )
            );
        }
        if (bf[0].order == bf[1].order) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument detected! You should not set the same order of boundary mode "
                    "twice: bf1.order == {0:d} while bf2.order == {1:d}",
                    bf[0].order, bf[1].order
                )
            );
        }
#endif
        constexpr std::array<param_type, 3> kFactors{-(7.0 / 60.0), -(2.0 / 15.0), -(7.0 / 60.0)};

        if (b0[0].order > b0[1].order) {
            std::swap(b0[0], b0[1]);
        }
        if (bf[0].order > bf[1].order) {
            std::swap(bf[0], bf[1]);
        }

        const std::size_t size{m_ts.size()};
        std::vector<param_type> hs(size - 1);
        std::vector<value_type> ds(size - 1);
        std::vector<param_type> a_diag(size * 2);
        std::vector<param_type> a_low_1(size * 2 - 1);
        std::vector<param_type> a_low_2(size + 1);
        std::vector<param_type> a_low_3(size);
        std::vector<param_type> a_low_4(size - 1);
        std::vector<param_type> a_up_1(size * 2 - 1);
        std::vector<param_type> a_up_2(size + 1);
        std::vector<param_type> a_up_3(size);
        std::vector<param_type> a_up_4(size - 1);
        std::vector<value_type> b(size * 2, value_type{0.0});

        hs[0] = m_ts[1] - m_ts[0];
        ds[0] = (m_ys[1] - m_ys[0]) / hs[0];
        a_low_2[0] = 0.0;
        a_up_2[0] = 0.0;
        for (std::size_t i{1}; i < size - 1; ++i) {
            hs[i] = m_ts[i + 1] - m_ts[i];
            ds[i] = (m_ys[i + 1] - m_ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            a_low_2[i + 1] = -(6.0 / hs[i]);
            a_low_3[i] = 6.0 / hs[i] + 6.0 / hs[i - 1];
            a_low_4[i - 1] = -(6.0 / hs[i - 1]);
            a_up_2[i] = hs[i - 1] * hs[i - 1] * hs[i - 1] * kFactors[0];
            a_up_3[i] = (hs[i - 1] * hs[i - 1] * hs[i - 1] + hs[i] * hs[i] * hs[i]) * kFactors[1];
            a_up_4[i] = hs[i] * hs[i] * hs[i] * kFactors[2];
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }
        a_low_2[size] = 0.0;
        a_up_2[size] = 0.0;

        std::copy(a_diag.cbegin(), a_diag.cbegin() + size, a_diag.begin() + size);
        std::copy(hs.cbegin(), hs.cend(), a_low_1.begin());
        std::copy(hs.cbegin(), hs.cend(), a_low_1.begin() + size);
        std::copy(hs.cbegin(), hs.cend(), a_up_1.begin());
        std::copy(hs.cbegin(), hs.cend(), a_up_1.begin() + size);

        a_low_1[size - 1] = 0.0;
        a_up_1[size - 1] = 0.0;

        if (b0[0].order == 1 && b0[1].order == 2) {
            const param_type h2 = hs[0] * hs[0];
            a_diag[size] = 1.0;
            a_low_2[1] = -7.5 / h2;
            a_low_3[0] = -15.0 / h2;
            a_up_1[size] = 0.875;
            b[size] = (b0[0].derivative - ds[0]) / (h2 * hs[0]) * 45.0;
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[1].derivative;
        } else if (b0[0].order == 2 && b0[1].order == 4) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[0].derivative;
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        } else if (b0[0].order == 1 && b0[1].order == 3) {
            const param_type h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[0].derivative) / hs[0] * 3.0;
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / h2);
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[1].derivative / hs[0] * 3.0);
        } else if (b0[0].order == 2 && b0[1].order == 3) {
            a_diag[0] = 1.0;
            a_up_1[0] = 0.0;
            a_up_3[0] = 0.0;
            a_up_4[0] = 0.0;
            b[0] = b0[0].derivative;
            a_diag[size] = 1.0;
            a_low_2[1] = -(3.0 / (hs[0] * hs[0]));
            a_low_3[0] = -a_low_2[1];
            a_up_1[size] = 0.5;
            b[size] = -(b0[1].derivative / hs[0] * 3.0);
        } else if (b0[0].order == 3 && b0[1].order == 4) {
            const param_type h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = -1.0;
            a_up_3[0] = h2 / 3.0;
            a_up_4[0] = h2 / 6.0;
            b[0] = -b0[0].derivative * hs[0];
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        } else {
            const param_type h2 = hs[0] * hs[0];
            a_diag[0] = 1.0;
            a_up_1[0] = 0.5;
            a_up_3[0] = h2 * kFactors[1] * 0.5;
            a_up_4[0] = h2 * kFactors[2] * 0.5;
            b[0] = (ds[0] - b0[0].derivative) / hs[0] * 3.0;
            a_diag[size] = 1.0;
            a_low_2[1] = 0.0;
            a_low_3[0] = 0.0;
            a_up_1[size] = 0.0;
            b[size] = b0[1].derivative;
        }

        if (bf[0].order == 1 && bf[1].order == 2) {
            const param_type h2 = hs[size - 2] * hs[size - 2];
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.875;
            a_low_3[size - 1] = -15.0 / h2;
            a_low_4[size - 2] = -7.5 / h2;
            b[size * 2 - 1] = (ds[size - 2] - bf[0].derivative) / (h2 * hs[size - 2]) * 45.0;
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[1].derivative;
        } else if (bf[0].order == 2 && bf[1].order == 4) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[0].derivative;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        } else if (bf[0].order == 1 && bf[1].order == 3) {
            const param_type h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[0].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / h2;
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[1].derivative / hs[size - 2] * 3.0;
        } else if (bf[0].order == 2 && bf[1].order == 3) {
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.0;
            a_up_2[size - 1] = 0.0;
            a_up_3[size - 1] = 0.0;
            b[size - 1] = bf[0].derivative;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.5;
            a_low_3[size - 1] = 3.0 / (hs[size - 2] * hs[size - 2]);
            a_low_4[size - 2] = -a_low_3[size - 1];
            b[size * 2 - 1] = bf[1].derivative / hs[size - 2] * 3.0;
        } else if (bf[0].order == 3 && bf[1].order == 4) {
            const param_type h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = -1.0;
            a_up_2[size - 1] = h2 / 6.0;
            a_up_3[size - 1] = h2 / 3.0;
            b[size - 1] = bf[0].derivative * hs[size - 2];
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        } else {
            const param_type h2 = hs[size - 2] * hs[size - 2];
            a_diag[size - 1] = 1.0;
            a_low_1[size - 2] = 0.5;
            a_up_2[size - 1] = h2 * kFactors[0] * 0.5;
            a_up_3[size - 1] = h2 * kFactors[1] * 0.5;
            b[size - 1] = (bf[0].derivative - ds[size - 2]) / hs[size - 2] * 3.0;
            a_diag[size * 2 - 1] = 1.0;
            a_low_1[size * 2 - 2] = 0.0;
            a_low_3[size - 1] = 0.0;
            a_low_4[size - 2] = 0.0;
            b[size * 2 - 1] = bf[1].derivative;
        }

        const OutriggerMatrix A{std::move(a_low_4), std::move(a_low_3), std::move(a_low_2),
                                std::move(a_low_1), std::move(a_diag),  std::move(a_up_1),
                                std::move(a_up_2),  std::move(a_up_3),  std::move(a_up_4)};

        const std::vector<value_type> x = A.gaussSeidel(b);

        m_ddys.resize(size);
        m_d4ys.resize(size);
        std::copy(x.cbegin(), x.cbegin() + size, m_ddys.begin());
        std::copy(x.cbegin() + size, x.cend(), m_d4ys.begin());
    }

    explicit PiecewiseQuinticFunction1(
        [[maybe_unused]] periodic_tag tag, std::vector<param_type> ts, std::vector<value_type> ys
    )
        : m_ts{std::move(ts)}, m_ys{std::move(ys)} {
#if BOYLE_CHECK_PARAMS == 1
        if (m_ts.size() < 2 || m_ys.size() < 2) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! sizes of ts, ys must be greater than 2: ts.size() "
                    "= "
                    "{0:d} while ys.size() = {1:d}",
                    m_ts.size(), m_ys.size()
                )
            );
        }
        if (m_ts.size() != m_ys.size()) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! ts, ys must share the same size: ts.size() = "
                    "{0:d} "
                    "while ys.size() = {1:d}",
                    m_ts.size(), m_ys.size()
                )
            );
        }
        if (!std::is_sorted(m_ts.cbegin(), m_ts.cend())) [[unlikely]] {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts has to be a sorted array!")
            );
        }
        if (hasDuplicates(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, kDuplicateCriterion))
            [[unlikely]] {
            throw std::invalid_argument(
                std::format("Invalid arguments detected! ts can not have duplicated elements!")
            );
        }
        if (std::abs(m_ys.front() - m_ys.back()) > kDuplicateCriterion) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid arguments detected! When choosing periodic boundary condition, it "
                    "requires ys.front() == m_ys.back() while m_ys.front() = {0}, m_ys.back() = "
                    "{1} "
                    "here.",
                    m_ys.front(), m_ys.back()
                )
            );
        }
#endif
        constexpr std::array<param_type, 3> kFactors{-(7.0 / 60.0), -(2.0 / 15.0), -(7.0 / 60.0)};

        const std::size_t size{m_ts.size() - 1};
        std::vector<param_type> hs(size);
        std::vector<value_type> ds(size);
        std::vector<param_type> a_diag(size * 2);
        std::vector<param_type> a_low_1(size * 2 - 1);
        std::vector<param_type> a_low_2(size + 1);
        std::vector<param_type> a_low_3(size);
        std::vector<param_type> a_low_4(size - 1);
        param_type a_bottom;
        std::vector<param_type> a_up_1(size * 2 - 1);
        std::vector<param_type> a_up_2(size + 1);
        std::vector<param_type> a_up_3(size);
        std::vector<param_type> a_up_4(size - 1);
        param_type a_top;
        std::vector<value_type> b(size * 2, value_type{0.0});

        hs[0] = m_ts[1] - m_ts[0];
        ds[0] = (m_ys[1] - m_ys[0]) / hs[0];
        for (std::size_t i{1}; i < size - 1; ++i) {
            hs[i] = m_ts[i + 1] - m_ts[i];
            ds[i] = (m_ys[i + 1] - m_ys[i]) / hs[i];
            a_diag[i] = (hs[i] + hs[i - 1]) * 2.0;
            a_low_2[i + 1] = -(6.0 / hs[i]);
            a_low_3[i] = 6.0 / hs[i] + 6.0 / hs[i - 1];
            a_low_4[i - 1] = -(6.0 / hs[i - 1]);
            a_up_2[i] = hs[i - 1] * hs[i - 1] * hs[i - 1] * kFactors[0];
            a_up_3[i] = (hs[i - 1] * hs[i - 1] * hs[i - 1] + hs[i] * hs[i] * hs[i]) * kFactors[1];
            a_up_4[i] = hs[i] * hs[i] * hs[i] * kFactors[2];
            b[i] = (ds[i] - ds[i - 1]) * 6.0;
        }
        hs[size - 1] = m_ts[size] - m_ts[size - 1];
        ds[size - 1] = (m_ys[size] - m_ys[size - 1]) / hs[size - 1];
        a_diag[size - 1] = (hs[size - 1] + hs[size - 2]) * 2.0;
        a_low_2[size] = -(6.0 / hs[size - 1]);
        a_low_3[size - 1] = 6.0 / hs[size - 1] + 6.0 / hs[size - 2];
        a_low_4[size - 2] = -(6.0 / hs[size - 2]);
        a_up_2[size - 1] = hs[size - 2] * hs[size - 2] * hs[size - 2] * kFactors[0];
        a_up_3[size - 1] = (hs[size - 2] * hs[size - 2] * hs[size - 2] +
                            hs[size - 1] * hs[size - 1] * hs[size - 1]) *
                           kFactors[1];
        b[size - 1] = (ds[size - 1] - ds[size - 2]) * 6.0;

        a_diag[0] = (hs[0] + hs[size - 1]) * 2.0;
        a_low_3[0] = 6.0 / hs[0] + 6.0 / hs[size - 1];
        a_up_3[0] =
            (hs[size - 1] * hs[size - 1] * hs[size - 1] + hs[0] * hs[0] * hs[0]) * kFactors[1];
        a_up_4[0] = hs[0] * hs[0] * hs[0] * kFactors[2];
        a_low_2[1] = -(6.0 / hs[0]);
        b[0] = (ds[0] - ds[size - 1]) * 6.0;

        a_low_2[0] = hs[size - 1];
        a_up_2[0] = hs[size - 1];
        a_low_2[size] = hs[size - 1];
        a_up_2[size] = hs[size - 1];
        a_bottom = -(6.0 / hs[size - 1]);
        a_low_1[size - 1] = -(6.0 / hs[size - 1]);
        a_up_1[size - 1] = hs[size - 1] * hs[size - 1] * hs[size - 1] * kFactors[2];
        a_top = hs[size - 1] * hs[size - 1] * hs[size - 1] * kFactors[0];

        std::copy(a_diag.cbegin(), a_diag.cbegin() + size, a_diag.begin() + size);
        std::copy(hs.cbegin(), hs.cend() - 1, a_low_1.begin());
        std::copy(hs.cbegin(), hs.cend() - 1, a_low_1.begin() + size);
        std::copy(hs.cbegin(), hs.cend() - 1, a_up_1.begin());
        std::copy(hs.cbegin(), hs.cend() - 1, a_up_1.begin() + size);

        const PeriodicOutriggerMatrix A{
            a_bottom,
            std::move(a_low_4),
            std::move(a_low_3),
            std::move(a_low_2),
            std::move(a_low_1),
            std::move(a_diag),
            std::move(a_up_1),
            std::move(a_up_2),
            std::move(a_up_3),
            std::move(a_up_4),
            a_top
        };

        const std::vector<value_type> x = A.gaussSeidel(b);

        m_ddys.resize(size + 1);
        m_d4ys.resize(size + 1);
        std::copy(x.cbegin(), x.cbegin() + size, m_ddys.begin());
        std::copy(x.cbegin() + size, x.cend(), m_d4ys.begin());
        m_ddys[size] = m_ddys[0];
        m_d4ys[size] = m_d4ys[0];
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto eval(param_type t) const noexcept -> value_type {
        constexpr std::array<param_type, 4> kFactors{
            -(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0
        };
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            const param_type h{m_ts[1] - m_ts[0]};
            const param_type ratio{(t - m_ts[0]) / h};
            return lerp(m_ys[0], m_ys[1], ratio) +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * (t - m_ts[0]) * h +
                   (m_d4ys[0] * kFactors[2] + m_d4ys[1] * kFactors[3]) * (t - m_ts[0]) * h * h * h;
        }
        if (pos == m_ts.size()) {
            const param_type h{m_ts[pos - 2] - m_ts[pos - 1]};
            const param_type ratio{(t - m_ts[pos - 1]) / h};
            return lerp(m_ys[pos - 1], m_ys[pos - 2], ratio) +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) *
                       (t - m_ts[pos - 1]) * h +
                   (m_d4ys[pos - 1] * kFactors[2] + m_d4ys[pos - 2] * kFactors[3]) *
                       (t - m_ts[pos - 1]) * h * h * h;
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative(param_type t) const noexcept -> value_type {
        constexpr std::array<param_type, 4> kFactors{
            -(1.0 / 3.0), -(1.0 / 6.0), 1.0 / 45.0, 7.0 / 360.0
        };
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            const param_type h{m_ts[1] - m_ts[0]};
            return (m_ys[1] - m_ys[0]) / h +
                   (m_ddys[0] * kFactors[0] + m_ddys[1] * kFactors[1]) * h +
                   (m_d4ys[0] * kFactors[2] + m_d4ys[1] * kFactors[3]) * h * h * h;
        }
        if (pos == m_ts.size()) {
            const param_type h{m_ts[pos - 2] - m_ts[pos - 1]};
            return (m_ys[pos - 2] - m_ys[pos - 1]) / h +
                   (m_ddys[pos - 1] * kFactors[0] + m_ddys[pos - 2] * kFactors[1]) * h +
                   (m_d4ys[pos - 1] * kFactors[2] + m_d4ys[pos - 2] * kFactors[3]) * h * h * h;
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio{(t - m_ts[pos - 1]) / h};
        return quinerpd(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure, always_inline]]
    auto derivative(param_type t, unsigned int order) const noexcept(!BOYLE_CHECK_PARAMS)
        -> value_type {
#if BOYLE_CHECK_PARAMS == 1
        if (order < 1 || order > 5) [[unlikely]] {
            throw std::invalid_argument(
                std::format(
                    "Invalid argument error! The PiecewiseLinearFunction only has 1st, 2nd, 3rd, "
                    "4th "
                    "order derivatives: order = {0:d}.",
                    order
                )
            );
        }
#endif
        value_type result;
        switch (order) {
        case 1:
            result = derivative(t);
            break;
        case 2:
            result = derivative2(t);
            break;
        case 3:
            result = derivative3(t);
            break;
        case 4:
            result = derivative4(t);
            break;
        case 5:
            result = derivative5(t);
            break;
        default:
            std::unreachable();
        }
        return result;
    }

    [[using gnu: pure]]
    auto integral(param_type lower_bound, param_type upper_bound) const noexcept -> value_type {
        constexpr std::array<param_type, 2> kFactors{-(1.0 / 24.0), 1.0 / 240.0};
        int sign{1};
        if (lower_bound > upper_bound) {
            std::swap(lower_bound, upper_bound);
            sign = -1;
        }
        const std::size_t size{m_ts.size()};
        const std::size_t istart =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), lower_bound) - m_ts.cbegin();
        const std::size_t iend =
            std::lower_bound(m_ts.cbegin(), m_ts.cend(), upper_bound) - m_ts.cbegin();
        param_type h, h3;
        if (istart == size || iend == 0 || istart == iend) {
            h = upper_bound - lower_bound;
            return (eval(lower_bound) + eval(upper_bound)) * h * 0.5 * sign;
        }
        value_type result{0.0};
        h = m_ts[istart] - lower_bound;
        h3 = h * h * h;
        result += (eval(lower_bound) + m_ys[istart]) * h * 0.5 +
                  (derivative2(lower_bound) + m_ddys[istart]) * h3 * kFactors[0] +
                  (derivative4(lower_bound) + m_d4ys[istart]) * h3 * h * h * kFactors[1];
        for (std::size_t i{istart}; i < iend - 1; ++i) {
            h = m_ts[i + 1] - m_ts[i];
            h3 = h * h * h;
            result += (m_ys[i] + m_ys[i + 1]) * h * 0.5 +
                      (m_ddys[i] + m_ddys[i + 1]) * h3 * kFactors[0] +
                      (m_d4ys[i] + m_d4ys[i + 1]) * h3 * h * h * kFactors[1];
        }
        h = upper_bound - m_ts[iend - 1];
        h3 = h * h * h;
        result += (m_ys[iend - 1] + eval(upper_bound)) * h * 0.5 +
                  (m_ddys[iend - 1] + derivative2(upper_bound)) * h3 * kFactors[0] +
                  (m_d4ys[iend - 1] + derivative4(upper_bound)) * h3 * h * h * kFactors[1];
        result *= sign;
        return result;
    }

    [[using gnu: pure, always_inline]]
    auto minT() const noexcept -> param_type {
        return m_ts.front();
    }

    [[using gnu: pure, always_inline]]
    auto maxT() const noexcept -> param_type {
        return m_ts.back();
    }

    [[using gnu: pure]]
    auto minY() const noexcept -> T
        requires std::floating_point<value_type>
    {
        std::size_t pos = std::min_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        param_type h, ratio;
        value_type derivative, derivative2;
        if (pos == 0) {
            h = m_ts[1] - m_ts[0];
            derivative =
                quinerpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], m_d4ys[0], m_d4ys[1], 0.0, h);
            if (derivative < 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], 1.0, h
            );
            if (derivative < 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = quinerpd(
                m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], m_d4ys[pos],
                m_d4ys[pos + 1], 0.0, h
            );
            if (derivative < 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], ratio, h
            );
            derivative2 =
                cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
            ratio -= derivative / (derivative2 * h);
        }
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure]]
    auto maxY() const noexcept -> T
        requires std::floating_point<value_type>
    {
        std::size_t pos = std::max_element(m_ys.cbegin(), m_ys.cend()) - m_ys.cbegin();
        param_type h, ratio;
        value_type derivative, derivative2;
        if (pos == 0) {
            h = m_ts[1] - m_ts[0];
            derivative =
                quinerpd(m_ys[0], m_ys[1], m_ddys[0], m_ddys[1], m_d4ys[0], m_d4ys[1], 0.0, h);
            if (derivative > 0.0) {
                pos += 1;
            } else {
                return m_ys[0];
            }
        } else if (pos == m_ts.size() - 1) {
            h = m_ts[pos] - m_ts[pos - 1];
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], 1.0, h
            );
            if (derivative > 0.0) {
                return m_ys[pos];
            }
        } else {
            h = m_ts[pos + 1] - m_ts[pos];
            derivative = quinerpd(
                m_ys[pos], m_ys[pos + 1], m_ddys[pos], m_ddys[pos + 1], m_d4ys[pos],
                m_d4ys[pos + 1], 0.0, h
            );
            if (derivative > 0.0) {
                pos += 1;
            }
        }
        h = m_ts[pos] - m_ts[pos - 1];
        ratio = 0.5;
        for (std::size_t num_iter{3}; num_iter != 0U; --num_iter) {
            derivative = quinerpd(
                m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1],
                m_d4ys[pos], ratio, h
            );
            derivative2 =
                cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
            ratio -= derivative / (derivative2 * h);
        }
        return quinerp(
            m_ys[pos - 1], m_ys[pos], m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos],
            ratio, h
        );
    }

    [[using gnu: pure, always_inline, hot]]
    auto operator()(param_type t) const noexcept -> value_type {
        return eval(t);
    }

    [[using gnu: pure, always_inline]]
    auto ts() const noexcept -> std::span<const param_type> {
        return m_ts;
    }

    [[using gnu: pure, always_inline]]
    auto ys() const noexcept -> std::span<const value_type> {
        return m_ys;
    }

    [[using gnu: pure, always_inline]]
    auto ddys() const noexcept -> std::span<const value_type> {
        return m_ddys;
    }

    [[using gnu: pure, always_inline]]
    auto d4ys() const noexcept -> std::span<const value_type> {
        return m_d4ys;
    }

  private:
    struct [[nodiscard]] OutriggerMatrix final {
        [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
        auto gaussSeidel(std::span<const T> b) const noexcept -> std::vector<value_type> {
            const std::size_t mat_size{a_diag.size()};
            const std::size_t half_mat_size{mat_size / 2};
            std::vector<value_type> x(mat_size, 0.0);
            for (std::size_t num_iter{40}; num_iter != 0U; --num_iter) {
                x[0] = (b[0] - a_up_1[0] * x[1] - a_up_2[0] * x[half_mat_size - 1] -
                        a_up_3[0] * x[half_mat_size] - a_up_4[0] * x[half_mat_size + 1]) /
                       a_diag[0];
                for (std::size_t i{1}; i < half_mat_size - 1; ++i) {
                    x[i] =
                        (b[i] - a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1] -
                         a_up_2[i] * x[half_mat_size + i - 1] - a_up_3[i] * x[half_mat_size + i] -
                         a_up_4[i] * x[half_mat_size + i + 1]) /
                        a_diag[i];
                }
                x[half_mat_size - 1] = (b[half_mat_size - 1] - a_low_2[0] * x[0] -
                                        a_low_1[half_mat_size - 2] * x[half_mat_size - 2] -
                                        a_up_1[half_mat_size - 1] * x[half_mat_size] -
                                        a_up_2[half_mat_size - 1] * x[mat_size - 2] -
                                        a_up_3[half_mat_size - 1] * x[mat_size - 1]) /
                                       a_diag[half_mat_size - 1];
                x[half_mat_size] = (b[half_mat_size] - a_low_3[0] * x[0] - a_low_2[1] * x[1] -
                                    a_low_1[half_mat_size - 1] * x[half_mat_size - 1] -
                                    a_up_1[half_mat_size] * x[half_mat_size + 1] -
                                    a_up_2[half_mat_size] * x[mat_size - 1]) /
                                   a_diag[half_mat_size];
                for (std::size_t i{half_mat_size + 1}; i < mat_size - 1; ++i) {
                    x[i] = (b[i] - a_low_4[i - half_mat_size - 1] * x[i - half_mat_size - 1] -
                            a_low_3[i - half_mat_size] * x[i - half_mat_size] -
                            a_low_2[i - half_mat_size + 1] * x[i - half_mat_size + 1] -
                            a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1]) /
                           a_diag[i];
                }
                x[mat_size - 1] =
                    (b[mat_size - 1] - a_low_4[half_mat_size - 2] * x[half_mat_size - 2] -
                     a_low_3[half_mat_size - 1] * x[half_mat_size - 1] -
                     a_low_2[half_mat_size] * x[half_mat_size] -
                     a_low_1[mat_size - 2] * x[mat_size - 2]) /
                    a_diag[mat_size - 1];
            }
            return x;
        }

        std::vector<param_type> a_low_4;
        std::vector<param_type> a_low_3;
        std::vector<param_type> a_low_2;
        std::vector<param_type> a_low_1;
        std::vector<param_type> a_diag;
        std::vector<param_type> a_up_1;
        std::vector<param_type> a_up_2;
        std::vector<param_type> a_up_3;
        std::vector<param_type> a_up_4;
    };

    struct [[nodiscard]] PeriodicOutriggerMatrix final {
        [[using gnu: pure, flatten, leaf, hot]] [[nodiscard]]
        auto gaussSeidel(std::span<const T> b) const noexcept -> std::vector<value_type> {
            const std::size_t mat_size{a_diag.size()};
            const std::size_t half_mat_size{mat_size / 2};
            std::vector<value_type> x(mat_size, 0.0);
            for (std::size_t num_iter{40}; num_iter != 0U; --num_iter) {
                x[0] = (b[0] - a_up_1[0] * x[1] - a_up_2[0] * x[half_mat_size - 1] -
                        a_up_3[0] * x[half_mat_size] - a_up_4[0] * x[half_mat_size + 1] -
                        a_top * x[mat_size - 1]) /
                       a_diag[0];
                for (std::size_t i{1}; i < half_mat_size - 1; ++i) {
                    x[i] =
                        (b[i] - a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1] -
                         a_up_2[i] * x[half_mat_size + i - 1] - a_up_3[i] * x[half_mat_size + i] -
                         a_up_4[i] * x[half_mat_size + i + 1]) /
                        a_diag[i];
                }
                x[half_mat_size - 1] = (b[half_mat_size - 1] - a_low_2[0] * x[0] -
                                        a_low_1[half_mat_size - 2] * x[half_mat_size - 2] -
                                        a_up_1[half_mat_size - 1] * x[half_mat_size] -
                                        a_up_2[half_mat_size - 1] * x[mat_size - 2] -
                                        a_up_3[half_mat_size - 1] * x[mat_size - 1]) /
                                       a_diag[half_mat_size - 1];
                x[half_mat_size] = (b[half_mat_size] - a_low_3[0] * x[0] - a_low_2[1] * x[1] -
                                    a_low_1[half_mat_size - 1] * x[half_mat_size - 1] -
                                    a_up_1[half_mat_size] * x[half_mat_size + 1] -
                                    a_up_2[half_mat_size] * x[mat_size - 1]) /
                                   a_diag[half_mat_size];
                for (std::size_t i{half_mat_size + 1}; i < mat_size - 1; ++i) {
                    x[i] = (b[i] - a_low_4[i - half_mat_size - 1] * x[i - half_mat_size - 1] -
                            a_low_3[i - half_mat_size] * x[i - half_mat_size] -
                            a_low_2[i - half_mat_size + 1] * x[i - half_mat_size + 1] -
                            a_low_1[i - 1] * x[i - 1] - a_up_1[i] * x[i + 1]) /
                           a_diag[i];
                }
                x[mat_size - 1] = (b[mat_size - 1] - a_bottom * x[0] -
                                   a_low_4[half_mat_size - 2] * x[half_mat_size - 2] -
                                   a_low_3[half_mat_size - 1] * x[half_mat_size - 1] -
                                   a_low_2[half_mat_size] * x[half_mat_size] -
                                   a_low_1[mat_size - 2] * x[mat_size - 2]) /
                                  a_diag[mat_size - 1];
            }
            return x;
        }

        param_type a_bottom;
        std::vector<param_type> a_low_4;
        std::vector<param_type> a_low_3;
        std::vector<param_type> a_low_2;
        std::vector<param_type> a_low_1;
        std::vector<param_type> a_diag;
        std::vector<param_type> a_up_1;
        std::vector<param_type> a_up_2;
        std::vector<param_type> a_up_3;
        std::vector<param_type> a_up_4;
        param_type a_top;
    };

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative2(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio = (t - m_ts[pos - 1]) / h;
        return cuberp(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative3(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        const param_type h{m_ts[pos] - m_ts[pos - 1]};
        const param_type ratio{(t - m_ts[pos - 1]) / h};
        return cuberpd(m_ddys[pos - 1], m_ddys[pos], m_d4ys[pos - 1], m_d4ys[pos], ratio, h);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative4(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        const param_type ratio{(t - m_ts[pos - 1]) / (m_ts[pos] - m_ts[pos - 1])};
        return lerp(m_d4ys[pos - 1], m_d4ys[pos], ratio);
    }

    [[using gnu: pure, flatten, leaf, hot]]
    auto derivative5(param_type t) const noexcept -> value_type {
        const std::size_t pos =
            nearestUpperElement(std::ranges::subrange{m_ts.cbegin(), m_ts.cend()}, t) -
            m_ts.cbegin();
        if (pos == 0) {
            return static_cast<value_type>(0.0);
        }
        if (pos == m_ts.size()) {
            return static_cast<value_type>(0.0);
        }
        return (m_d4ys[pos] - m_d4ys[pos - 1]) / (m_ts[pos] - m_ts[pos - 1]);
    }

    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_ts;
        archive & m_ys;
        archive & m_ddys;
        archive & m_d4ys;
        return;
    }

    std::vector<param_type> m_ts{};
    std::vector<value_type> m_ys{};
    std::vector<value_type> m_ddys{};
    std::vector<value_type> m_d4ys{};
};

using PiecewiseQuinticFunction1s = PiecewiseQuinticFunction1<float>;
using PiecewiseQuinticFunction1d = PiecewiseQuinticFunction1<double>;

} // namespace boyle::math
