/**
 * @file piecewise_cubic_curve3_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-12-24
 *
 * @copyright Copyright (c) 2024 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "boyle/math/curves/piecewise_cubic_curve.hpp"

#include <cmath>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#include "boyle/math/dense/vec3.hpp"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

namespace {

bool plot_graph{false};

auto createFigureHandle() noexcept -> matplot::figure_handle {
    matplot::figure_handle fig = matplot::figure();
    fig->size(1600, 1000);
    return fig;
}

} // namespace

namespace boyle::math {

TEST_CASE("LeftHelixTest") {
    const auto exact_helix = [](double theta) noexcept -> Vec3d {
        return Vec3d{2.0 * std::cos(theta), 2.0 * std::sin(theta), theta};
    };

    constexpr std::size_t kNumAnchors{81};
    constexpr double kStart{0.0};
    constexpr double kEnd{M_PI * 4.0};
    constexpr double kStep{(kEnd - kStart) / (kNumAnchors - 1)};

    std::vector<Vec3d> anchor_points(kNumAnchors);
    std::ranges::generate(
        anchor_points, [t = kStart - kStep, h = kStep]() mutable noexcept -> Vec3d {
            t += h;
            return Vec3d{2.0 * std::cos(t), 2.0 * std::sin(t), t};
        }
    );

    PiecewiseCubicCurve3d::BoundaryMode b0, bf;
    SUBCASE("Clamped") {
        b0 = PiecewiseCubicCurve3d::BoundaryMode{
            1, Vec3d{0.0, 2.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0)}
        };
        bf = PiecewiseCubicCurve3d::BoundaryMode{
            1, Vec3d{0.0, 2.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0)}
        };
    }
    SUBCASE("Natural") {
        b0 = PiecewiseCubicCurve3d::BoundaryMode{2, Vec3d{-0.4, 0.0, 0.0}};
        bf = PiecewiseCubicCurve3d::BoundaryMode{2, Vec3d{-0.4, 0.0, 0.0}};
    }

    const PiecewiseCubicCurve3d helix{anchor_points, b0, bf};
    CHECK_EQ(helix.minS(), 0.0);
    CHECK_EQ(helix.maxS(), doctest::Approx(4.0 * M_PI * std::sqrt(5.0)).epsilon(1E-6));

    for (std::size_t i{1}; i < kNumAnchors; ++i) {
        const double theta{kStart + (i - 0.5) * kStep};
        const double s{theta * std::sqrt(5.0)};
        const Vec3d point = exact_helix(theta);
        const SlvTripletd slv = helix.inverse(point);
        const Vec3d tangent{
            -std::sin(theta) * 2.0 / std::sqrt(5.0), std::cos(theta) * 2.0 / std::sqrt(5.0),
            1 / std::sqrt(5.0)
        };
        const Vec3d normal{-std::cos(theta), -std::sin(theta), 0.0};
        const Vec3d binormal{
            std::sin(theta) / std::sqrt(5.0), -std::cos(theta) / std::sqrt(5.0), 2 / std::sqrt(5.0)
        };
        CHECK_EQ(slv.s, doctest::Approx(s).epsilon(1E-6));
        CHECK_EQ(slv.l, doctest::Approx(0.0).epsilon(1E-5));
        CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(2E-7));
        CHECK(helix.tangent(s).identicalTo(tangent, 1E-4));
        CHECK(helix.normal(s).identicalTo(normal, 1E-4));
        CHECK(helix.binormal(s).identicalTo(binormal, 1E-4));
        CHECK_EQ(helix.curvature(s), doctest::Approx(2.0 / 5.0).epsilon(6E-4));
        CHECK_EQ(helix.torsion(s), doctest::Approx(1.0 / 5.0).epsilon(5E-4));
    }

    const std::array<Vec3d, 2> milestone_points{
        Vec3d{std::sqrt(3.0) * 0.5, 0.5, M_PI / 6.0},
        Vec3d{-1.5, std::sqrt(3.0) * 1.5, M_PI * 2.0 / 3.0}
    };

    SlvTripletd slv = helix.inverse(milestone_points[0]);
    CHECK_EQ(slv.s, doctest::Approx(M_PI / 6.0 * std::sqrt(5.0)).epsilon(4E-5));
    CHECK_EQ(slv.l, doctest::Approx(1.0).epsilon(3E-6));
    CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(2E-5));
    CHECK(helix(slv).identicalTo(milestone_points[0], 3E-7));

    slv = helix.inverse(milestone_points[1]);
    CHECK_EQ(slv.s, doctest::Approx(M_PI * 2.0 / 3.0 * std::sqrt(5.0)).epsilon(4E-5));
    CHECK_EQ(slv.l, doctest::Approx(-1.0).epsilon(2E-6));
    CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(3E-5));
    CHECK(helix(slv).identicalTo(milestone_points[1], 1E-7));

    if (plot_graph) {
        using namespace matplot;

        std::array<double, kNumAnchors> anchor_xs, anchor_ys, anchor_zs;
        std::ranges::generate(anchor_xs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].x;
        });
        std::ranges::generate(anchor_ys, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].y;
        });
        std::ranges::generate(anchor_zs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].z;
        });

        constexpr std::size_t kNumExact{401};
        std::array<double, kNumExact> exact_xs, exact_ys, exact_zs;
        std::ranges::generate(
            exact_xs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::cos(t);
            }
        );
        std::ranges::generate(
            exact_ys,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::sin(t);
            }
        );
        std::ranges::generate(
            exact_zs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return t;
            }
        );

        constexpr std::size_t kNumTest{80};
        std::array<double, kNumTest> test_xs, test_ys, test_zs;
        {
            const double start_s{helix.minS()}, end_s{helix.maxS()};
            const double step_s{(end_s - start_s) / kNumTest};
            for (std::size_t i{0}; i < kNumTest; ++i) {
                const Vec3d point = helix(start_s + step_s * (i + 0.5));
                test_xs[i] = point.x;
                test_ys[i] = point.y;
                test_zs[i] = point.z;
            }
        }

        constexpr std::size_t kNumIntpl{401};
        std::array<double, kNumIntpl> intpl_xs, intpl_ys, intpl_zs;
        {
            const double start_s{helix.minS()}, end_s{helix.maxS()};
            const double step_s{(end_s - start_s) / (kNumIntpl - 1)};
            for (std::size_t i{0}; i < kNumIntpl; ++i) {
                const Vec3d point = helix(start_s + step_s * i);
                intpl_xs[i] = point.x;
                intpl_ys[i] = point.y;
                intpl_zs[i] = point.z;
            }
        }

        auto fig = createFigureHandle();
        fig->title("LeftHelix");
        hold(on);
        grid(on);
        axis(equal);
        plot3(anchor_xs, anchor_ys, anchor_zs, "b.");
        plot3(exact_xs, exact_ys, exact_zs, "b:");
        plot3(test_xs, test_ys, test_zs, "r*");
        plot3(intpl_xs, intpl_ys, intpl_zs, "r:");
        fig->show();
    }
}

TEST_CASE("RightHelixTest") {
    const auto exact_helix = [](double theta) noexcept -> Vec3d {
        return Vec3d{2.0 * std::cos(-theta), 2.0 * std::sin(-theta), theta};
    };

    constexpr std::size_t kNumAnchors{81};
    constexpr double kStart{0.0};
    constexpr double kEnd{M_PI * 4.0};
    constexpr double kStep{(kEnd - kStart) / (kNumAnchors - 1)};

    std::vector<Vec3d> anchor_points(kNumAnchors);
    std::ranges::generate(
        anchor_points, [t = kStart - kStep, h = kStep]() mutable noexcept -> Vec3d {
            t += h;
            return Vec3d{2.0 * std::cos(-t), 2.0 * std::sin(-t), t};
        }
    );

    PiecewiseCubicCurve3d::BoundaryMode b0, bf;
    SUBCASE("Clamped") {
        b0 = PiecewiseCubicCurve3d::BoundaryMode{
            1, Vec3d{0.0, -2.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0)}
        };
        bf = PiecewiseCubicCurve3d::BoundaryMode{
            1, Vec3d{0.0, -2.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0)}
        };
    }
    SUBCASE("Natural") {
        b0 = PiecewiseCubicCurve3d::BoundaryMode{2, Vec3d{-0.4, 0.0, 0.0}};
        bf = PiecewiseCubicCurve3d::BoundaryMode{2, Vec3d{-0.4, 0.0, 0.0}};
    }

    const PiecewiseCubicCurve3d helix{anchor_points, b0, bf};
    CHECK_EQ(helix.minS(), 0.0);
    CHECK_EQ(helix.maxS(), doctest::Approx(4.0 * M_PI * std::sqrt(5.0)).epsilon(1E-6));

    for (std::size_t i{1}; i < kNumAnchors; ++i) {
        const double theta{kStart + (i - 0.5) * kStep};
        const double s{theta * std::sqrt(5.0)};
        const Vec3d point = exact_helix(theta);
        const SlvTripletd slv = helix.inverse(point);
        const Vec3d tangent{
            std::sin(-theta) * 2.0 / std::sqrt(5.0), -std::cos(-theta) * 2.0 / std::sqrt(5.0),
            1 / std::sqrt(5.0)
        };
        const Vec3d normal{-std::cos(-theta), -std::sin(-theta), 0.0};
        const Vec3d binormal{
            std::sin(-theta) / std::sqrt(5.0), -std::cos(-theta) / std::sqrt(5.0),
            -2 / std::sqrt(5.0)
        };
        CHECK_EQ(slv.s, doctest::Approx(s).epsilon(1E-6));
        CHECK_EQ(slv.l, doctest::Approx(0.0).epsilon(1E-5));
        CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(2E-7));
        CHECK(helix.tangent(s).identicalTo(tangent, 1E-4));
        CHECK(helix.normal(s).identicalTo(normal, 1E-4));
        CHECK(helix.binormal(s).identicalTo(binormal, 1E-4));
        CHECK_EQ(helix.curvature(s), doctest::Approx(2.0 / 5.0).epsilon(6E-4));
        CHECK_EQ(helix.torsion(s), doctest::Approx(-1.0 / 5.0).epsilon(5E-4));
    }

    const std::array<Vec3d, 2> milestone_points{
        Vec3d{std::sqrt(3.0) * 0.5, 0.5, M_PI * 11.0 / 6.0},
        Vec3d{-1.5, std::sqrt(3.0) * 1.5, M_PI * 4.0 / 3.0}
    };

    SlvTripletd slv = helix.inverse(milestone_points[0]);
    CHECK_EQ(slv.s, doctest::Approx(M_PI * 11.0 / 6.0 * std::sqrt(5.0)).epsilon(4E-5));
    CHECK_EQ(slv.l, doctest::Approx(1.0).epsilon(1E-5));
    CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(2E-3));
    CHECK(helix(slv).identicalTo(milestone_points[0], 1E-3));

    slv = helix.inverse(milestone_points[1]);
    CHECK_EQ(slv.s, doctest::Approx(M_PI * 4.0 / 3.0 * std::sqrt(5.0)).epsilon(4E-5));
    CHECK_EQ(slv.l, doctest::Approx(-1.0).epsilon(2E-6));
    CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(3E-5));
    CHECK(helix(slv).identicalTo(milestone_points[1], 1E-7));

    if (plot_graph) {
        using namespace matplot;

        std::array<double, kNumAnchors> anchor_xs, anchor_ys, anchor_zs;
        std::ranges::generate(anchor_xs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].x;
        });
        std::ranges::generate(anchor_ys, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].y;
        });
        std::ranges::generate(anchor_zs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].z;
        });

        constexpr std::size_t kNumExact{401};
        std::array<double, kNumExact> exact_xs, exact_ys, exact_zs;
        std::ranges::generate(
            exact_xs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::cos(-t);
            }
        );
        std::ranges::generate(
            exact_ys,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::sin(-t);
            }
        );
        std::ranges::generate(
            exact_zs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return t;
            }
        );

        constexpr std::size_t kNumTest{80};
        std::array<double, kNumTest> test_xs, test_ys, test_zs;
        {
            const double start_s{helix.minS()}, end_s{helix.maxS()};
            const double step_s{(end_s - start_s) / kNumTest};
            for (std::size_t i{0}; i < kNumTest; ++i) {
                const Vec3d point = helix(start_s + step_s * (i + 0.5));
                test_xs[i] = point.x;
                test_ys[i] = point.y;
                test_zs[i] = point.z;
            }
        }

        constexpr std::size_t kNumIntpl{401};
        std::array<double, kNumIntpl> intpl_xs, intpl_ys, intpl_zs;
        {
            const double start_s{helix.minS()}, end_s{helix.maxS()};
            const double step_s{(end_s - start_s) / (kNumIntpl - 1)};
            for (std::size_t i{0}; i < kNumIntpl; ++i) {
                const Vec3d point = helix(start_s + step_s * i);
                intpl_xs[i] = point.x;
                intpl_ys[i] = point.y;
                intpl_zs[i] = point.z;
            }
        }

        auto fig = createFigureHandle();
        fig->title("RightHelix");
        hold(on);
        grid(on);
        axis(equal);
        plot3(anchor_xs, anchor_ys, anchor_zs, "b.");
        plot3(exact_xs, exact_ys, exact_zs, "b:");
        plot3(test_xs, test_ys, test_zs, "r*");
        plot3(intpl_xs, intpl_ys, intpl_zs, "r:");
        fig->show();
    }
}

TEST_CASE("LissajousCurveTest") {
    const auto exact_lissajous_curve = [](double theta) noexcept -> Vec3d {
        return Vec3d{
            2.0 * std::sin(theta), 2.0 * std::sin(3.0 * theta + M_PI_4),
            2.0 * std::sin(5.0 * theta + M_PI_2)
        };
    };

    constexpr std::size_t kNumAnchors{101};
    constexpr double kStart{0.0};
    constexpr double kEnd{2.0 * M_PI};
    constexpr double kStep{(kEnd - kStart) / (kNumAnchors - 1)};

    std::vector<Vec3d> anchor_points(kNumAnchors);
    std::ranges::generate(
        anchor_points, [t = kStart - kStep, h = kStep]() mutable noexcept -> Vec3d {
            t += h;
            return Vec3d{
                2.0 * std::sin(t), 2.0 * std::sin(3.0 * t + M_PI_4),
                2.0 * std::sin(5.0 * t + M_PI_2)
            };
        }
    );

    const PiecewiseCubicCurve3d lissajous_curve{periodic_tag{}, anchor_points};

    const auto& arc_lengths{lissajous_curve.arcLengths()};

    for (std::size_t i{1}; i < kNumAnchors; ++i) {
        const Vec3d point = exact_lissajous_curve(kStart + (i - 0.5) * kStep);
        const SlvTripletd slv = lissajous_curve.inverse(point, arc_lengths[i - 1], arc_lengths[i]);
        CHECK_EQ(slv.l, doctest::Approx(0.0).epsilon(3E-2));
        CHECK_EQ(slv.v, doctest::Approx(0.0).epsilon(7E-3));
    }

    if (plot_graph) {
        using namespace matplot;

        std::array<double, kNumAnchors> anchor_xs, anchor_ys, anchor_zs;
        std::ranges::generate(anchor_xs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].x;
        });
        std::ranges::generate(anchor_ys, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].y;
        });
        std::ranges::generate(anchor_zs, [&anchor_points, i = -1]() mutable noexcept -> double {
            ++i;
            return anchor_points[i].z;
        });

        constexpr std::size_t kNumExact{1001};
        std::array<double, kNumExact> exact_xs, exact_ys, exact_zs;
        std::ranges::generate(
            exact_xs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::sin(t);
            }
        );
        std::ranges::generate(
            exact_ys,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::sin(3.0 * t + M_PI_4);
            }
        );
        std::ranges::generate(
            exact_zs,
            [t = kStart - (kEnd - kStart) / (kNumExact - 1),
             h = (kEnd - kStart) / (kNumExact - 1)]() mutable noexcept -> double {
                t += h;
                return 2.0 * std::sin(5.0 * t + M_PI_2);
            }
        );

        constexpr std::size_t kNumTest{100};
        std::array<double, kNumTest> test_xs, test_ys, test_zs;
        {
            const double start_s{lissajous_curve.minS()}, end_s{lissajous_curve.maxS()};
            const double step_s{(end_s - start_s) / kNumTest};
            for (std::size_t i{0}; i < kNumTest; ++i) {
                const Vec3d point = lissajous_curve(start_s + step_s * (i + 0.5));
                test_xs[i] = point.x;
                test_ys[i] = point.y;
                test_zs[i] = point.z;
            }
        }

        constexpr std::size_t kNumIntpl{1001};
        std::array<double, kNumIntpl> intpl_xs, intpl_ys, intpl_zs;
        {
            const double start_s{lissajous_curve.minS()}, end_s{lissajous_curve.maxS()};
            const double step_s{(end_s - start_s) / (kNumIntpl - 1)};
            for (std::size_t i{0}; i < kNumIntpl; ++i) {
                const Vec3d point = lissajous_curve(start_s + step_s * i);
                intpl_xs[i] = point.x;
                intpl_ys[i] = point.y;
                intpl_zs[i] = point.z;
            }
        }

        auto fig = createFigureHandle();
        fig->title("Lissajous Curve");
        hold(on);
        grid(on);
        axis(equal);
        plot3(anchor_xs, anchor_ys, anchor_zs, "b.");
        plot3(exact_xs, exact_ys, exact_zs, "b:");
        plot3(test_xs, test_ys, test_zs, "r*");
        plot3(intpl_xs, intpl_ys, intpl_zs, "r:");
        fig->show();
    }
}

TEST_CASE("Serialization") {
    [[maybe_unused]] const auto exact_lissajous_curve = [](double theta) noexcept -> Vec3d {
        return Vec3d{
            2.0 * std::sin(theta), 2.0 * std::sin(3.0 * theta + M_PI_4),
            2.0 * std::sin(5.0 * theta + M_PI_2)
        };
    };

    constexpr std::size_t kNumAnchors{101};
    constexpr double kStart{0.0};
    constexpr double kEnd{2.0 * M_PI};
    constexpr double kStep{(kEnd - kStart) / (kNumAnchors - 1)};

    std::vector<Vec3d> anchor_points(kNumAnchors);
    std::ranges::generate(
        anchor_points, [t = kStart - kStep, h = kStep]() mutable noexcept -> Vec3d {
            t += h;
            return Vec3d{
                2.0 * std::sin(t), 2.0 * std::sin(3.0 * t + M_PI_4),
                2.0 * std::sin(5.0 * t + M_PI_2)
            };
        }
    );

    const PiecewiseCubicCurve3d lissajous_curve{periodic_tag{}, anchor_points};

    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << lissajous_curve;

    PiecewiseCubicCurve3d other_lissajous_curve;

    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_lissajous_curve;

    const auto& arc_length{lissajous_curve.arcLengths()};
    const auto& other_arc_length{other_lissajous_curve.arcLengths()};
    const auto& other_anchor_points{other_lissajous_curve.anchorPoints()};

    for (std::size_t i = 0; i < kNumAnchors; ++i) {
        CHECK_EQ(other_arc_length[i], doctest::Approx(arc_length[i]).epsilon(kEpsilon));
        CHECK_EQ(other_anchor_points[i].x, doctest::Approx(anchor_points[i].x).epsilon(kEpsilon));
        CHECK_EQ(other_anchor_points[i].y, doctest::Approx(anchor_points[i].y).epsilon(kEpsilon));
    }
}

} // namespace boyle::math

auto main(int argc, const char* argv[]) -> int {
    cxxopts::Options options(
        "piecewise_linear_curve3_test", "unit test of PiecewiseLinearCurve3 class"
    );
    options.add_options()(
        "plot-graph", "plot test graph", cxxopts::value<bool>()->default_value("false")
    );
    cxxopts::ParseResult result = options.parse(argc, argv);
    plot_graph = result["plot-graph"].as<bool>();
    doctest::Context context(argc, argv);
    return context.run();
}
