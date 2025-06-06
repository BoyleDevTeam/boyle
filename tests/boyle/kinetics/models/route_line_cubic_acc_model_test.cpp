/**
 * @file route_line_cubic_acc_model_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-18
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/kinetics/models/route_line_cubic_acc_model.hpp"

#include "cxxopts.hpp"
#include "matplot/matplot.h"

#include "boyle/math/utils.hpp"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

namespace {

bool plot_graph{false};

auto createFigureHandle() noexcept -> matplot::figure_handle {
    matplot::figure_handle fig = matplot::figure();
    fig->size(1600, 1000);
    return fig;
}

auto createAxesHandles(matplot::figure_handle& fig) noexcept -> std::vector<matplot::axes_handle> {
    std::vector<matplot::axes_handle> axes;
    axes.reserve(5);
    for (std::size_t i = 0; i < 5; ++i) {
        axes.emplace_back(fig->add_subplot(2, 3, i));
        axes[i]->grid(matplot::on);
        axes[i]->hold(matplot::on);
        axes[i]->font_size(12.0);
        axes[i]->x_axis().label_font_size(12.0);
        axes[i]->y_axis().label_font_size(12.0);
    }
    axes[0]->xlabel("t");
    axes[0]->ylabel("s");
    axes[1]->xlabel("t");
    axes[1]->ylabel("v");
    axes[2]->xlabel("t");
    axes[2]->ylabel("a");
    axes[3]->xlabel("t");
    axes[3]->ylabel("j");
    axes[4]->xlabel("t");
    axes[4]->ylabel("snap");
    return axes;
}

} // namespace

namespace boyle::kinetics {

TEST_CASE("CostTest") {
    auto s_profile = [](double t) -> double {
        return 14.0 + 5.0 * t + 2.8 * t * t + 0.5 * t * t * t;
    };
    auto v_profile = [](double t) -> double { return 5.0 + 5.6 * t + 1.5 * t * t; };

    constexpr std::size_t num_samples = 21;
    const std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{sample_ts};

    std::vector<double> state_vec(num_samples * 2);
    for (std::size_t i{0}; i < num_samples; ++i) {
        state_vec[i] = s_profile(sample_ts[i]);
        state_vec[num_samples + i] = v_profile(sample_ts[i]);
    }

    double cost{0.0};

    SUBCASE("VelocityCost") {
        route_line_acc_model.setVelocityCost(5.0, 20.0);
        cost = route_line_acc_model.qp_problem().cost({state_vec.cbegin(), state_vec.cend()});
        CHECK_EQ(
            cost,
            doctest::Approx(2195626.66667 * 20.0 / 20.0 - 500.0).epsilon(::boyle::math::kEpsilon)
        );
    }

    SUBCASE("AccelCost") {
        route_line_acc_model.setAccelCost(10.0);
        cost = route_line_acc_model.qp_problem().cost({state_vec.cbegin(), state_vec.cend()});
        CHECK_EQ(cost, doctest::Approx(31347.2 * 10.0 / 20.0).epsilon(math::kEpsilon));
    }

    SUBCASE("JerkCost") {
        route_line_acc_model.setJerkCost(100.0);
        cost = route_line_acc_model.qp_problem().cost({state_vec.cbegin(), state_vec.cend()});
        CHECK_EQ(cost, doctest::Approx(180.0 * 100.0 / 20.0).epsilon(1E-6));
    }

    SUBCASE("SoftFence") {
        std::vector<SoftFence1d> soft_fences{SoftFence1d(
            0, ::boyle::kinetics::Actio::PUSHING, {0.0, 20.0}, {-946.0, 3524.0}, 100.0, 10.0
        )};
        for (std::size_t i{0}; i < num_samples; ++i) {
            state_vec.push_back(std::max(223.5 * sample_ts[i] - 946.0 - state_vec[i], 0.0));
        }
        route_line_acc_model.setSoftFences(soft_fences);
        cost = route_line_acc_model.qp_problem().cost({state_vec.cbegin(), state_vec.cend()});
        CHECK_EQ(
            cost,
            doctest::Approx(2966.67 * 100.0 / 20.0 + 1059109.523809519 * 10.0 / 20.0).epsilon(1E-3)
        );
    }
}

TEST_CASE("TrivialScene") {
    const std::size_t num_samples = 21;
    std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);
    const Motion1d motion = route_line_acc_model.solve().first;

    CHECK_EQ(motion.s(0.0), doctest::Approx(0.0).epsilon(1E-2));
    CHECK_EQ(motion.s(10.0), doctest::Approx(50.0).epsilon(1E1));
    CHECK_EQ(motion.s(20.0), doctest::Approx(100.0).epsilon(1E-2));

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = createFigureHandle();
        std::vector<axes_handle> axes = createAxesHandles(fig);
        axes[0]->plot(plot_ts, plot_ss, ".");
        axes[1]->plot(plot_ts, plot_vs, ".");
        axes[2]->plot(plot_ts, plot_as, ".");
        axes[3]->plot(plot_ts, plot_js, ".");
        axes[4]->plot(plot_ts, plot_snaps, ".");
        fig->show();
    }
}

TEST_CASE("FenceScene") {
    const std::size_t num_samples = 21;
    std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    std::vector<double> bound_ts_1{8.1, 8.9};
    std::vector<double> bound_ss_1{60.0, 60.0};
    std::vector<double> bound_ts_2{3.2, 3.8};
    std::vector<double> bound_ss_2{20.0, 20.0};

    Motion1d motion;

    SUBCASE("SoftFences") {
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d{
                0, ::boyle::kinetics::Actio::BLOCKING, bound_ts_1, bound_ss_1, 1000.0, 10.0
            },
            SoftFence1d{1, ::boyle::kinetics::Actio::PUSHING, bound_ts_2, bound_ss_2, 1000.0, 10.0}
        };

        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().first;
    }

    SUBCASE("HardFences") {
        std::vector<HardFence1d> hard_fences{
            HardFence1d{0, ::boyle::kinetics::Actio::BLOCKING, bound_ts_1, bound_ss_1},
            HardFence1d{1, ::boyle::kinetics::Actio::PUSHING, bound_ts_2, bound_ss_2}
        };

        route_line_acc_model.setHardFences(hard_fences);
        motion = route_line_acc_model.solve().first;

        CHECK_LT(motion.s(bound_ts_1.back() - 1E-2), bound_ss_1.back());
        CHECK_GT(motion.s(bound_ts_2.front() + 1E-2), bound_ss_2.front());
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = createFigureHandle();
        std::vector<axes_handle> axes = createAxesHandles(fig);
        axes[0]->plot(plot_ts, plot_ss, ".");
        axes[0]->plot(bound_ts_1, bound_ss_1, "r-");
        axes[0]->plot({bound_ts_1.front(), bound_ts_1.front()}, {bound_ss_1.front(), 100.0}, "r-");
        axes[0]->plot({bound_ts_1.back(), bound_ts_1.back()}, {bound_ss_1.back(), 100.0}, "r-");
        axes[0]->plot(bound_ts_2, bound_ss_2, "r-");
        axes[0]->plot({bound_ts_2.front(), bound_ts_2.front()}, {0.0, bound_ss_2.front()}, "r-");
        axes[0]->plot({bound_ts_2.back(), bound_ts_2.back()}, {0.0, bound_ss_2.back()}, "r-");
        axes[1]->plot(plot_ts, plot_vs, ".");
        axes[2]->plot(plot_ts, plot_as, ".");
        axes[3]->plot(plot_ts, plot_js, ".");
        axes[4]->plot(plot_ts, plot_snaps, ".");
        show();
    }
}

TEST_CASE("HeadwayScene") {
    const std::size_t num_samples = 51;
    std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 10.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    Motion1d motion;
    std::vector<double> bound_ts;
    std::vector<double> bound_ss;

    SUBCASE("Yield") {
        route_line_acc_model.setInitialState(0.0, 20.0, 0.0);
        route_line_acc_model.setVelocityCost(20.0, 1.0);
        bound_ts = {0.0, 10.0};
        bound_ss = {20.0, 120.0};
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::kinetics::Actio::BLOCKING, bound_ts, bound_ss, 10.0, 10.0)
        };
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().first;
    }

    SUBCASE("Overtake") {
        route_line_acc_model.setInitialState(0.0, 10.0, 0.0);
        route_line_acc_model.setVelocityCost(10.0, 1.0);
        bound_ts = {0.0, 10.0};
        bound_ss = {-10.0, 140.0};
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::kinetics::Actio::PUSHING, bound_ts, bound_ss, 10.0, 10.0)
        };
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().first;
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 10.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = createFigureHandle();
        std::vector<axes_handle> axes = createAxesHandles(fig);
        axes[0]->plot({bound_ts}, {bound_ss}, "r-");
        axes[0]->plot(plot_ts, plot_ss, "-")->line_width(2.0);
        axes[1]->plot(plot_ts, plot_vs, "-")->line_width(2.0);
        axes[2]->plot(plot_ts, plot_as, "-")->line_width(2.0);
        axes[3]->plot(plot_ts, plot_js, "-")->line_width(2.0);
        axes[4]->plot(plot_ts, plot_snaps, "-")->line_width(2.0);
        show();
    }
}

} // namespace boyle::kinetics

auto main(int argc, const char* argv[]) -> int {
    cxxopts::Options options(
        "piecewise_linear_function_test", "unit test of PiecewiseLinearFunction class"
    );
    options.add_options()(
        "plot-graph", "plot test graph", cxxopts::value<bool>()->default_value("false")
    );
    cxxopts::ParseResult result = options.parse(argc, argv);
    plot_graph = result["plot-graph"].as<bool>();
    doctest::Context context(argc, argv);
    return context.run();
}
