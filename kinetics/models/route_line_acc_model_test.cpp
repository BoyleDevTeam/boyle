/**
 * @file route_line_acc_model_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-18
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#include "kinetics/models/route_line_acc_model.h"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "math/utils.hpp"

namespace {

bool plot_graph;

} // namespace

namespace tiny_pnc {
namespace kinetics {

TEST_CASE("TrivialScene") {
    constexpr std::size_t num_samples = 41;
    const std::vector<double> sample_ts = tiny_pnc::math::linspace(0.0, 20.0, num_samples);
    RouteLineAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(1.0);
    route_line_acc_model.setJerkCost(10.0);
    route_line_acc_model.setSnapCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);
    RouteLineAccModel::Result result = route_line_acc_model.solve();
    const Motion1d& motion = result.motion;

    CHECK_EQ(motion.s(0.0), doctest::Approx(0.0).epsilon(1E-2));
    CHECK_EQ(motion.s(10.0), doctest::Approx(50.0).epsilon(1E-2));
    CHECK_EQ(motion.s(20.0), doctest::Approx(100.0).epsilon(1E-2));

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = tiny_pnc::math::linspace(
            0.0 + tiny_pnc::math::kEpsilon, 20.0 - tiny_pnc::math::kEpsilon, num_samples
        );
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
        {
            std::vector<axes_handle> ax;
            for (std::size_t i = 0; i < 5; ++i) {
                ax.emplace_back(subplot(2, 3, i));
            }
            ax[0]->grid(on);
            ax[0]->plot(plot_ts, plot_ss, ".");

            ax[1]->grid(on);
            ax[1]->plot(plot_ts, plot_vs, ".");

            ax[2]->grid(on);
            ax[2]->plot(plot_ts, plot_as, ".");

            ax[3]->grid(on);
            ax[3]->plot(plot_ts, plot_js, ".");

            ax[4]->grid(on);
            ax[4]->plot(plot_ts, plot_snaps, ".");
            show();

            for (std::size_t i = 0; i < 5; ++i) {
                ax[i]->clear();
            }
        }
    }
}

TEST_CASE("FenceScene") {
    constexpr std::size_t num_samples = 41;
    const std::vector<double> sample_ts = tiny_pnc::math::linspace(0.0, 20.0, num_samples);
    RouteLineAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(1.0);
    route_line_acc_model.setJerkCost(10.0);
    route_line_acc_model.setSnapCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    tiny_pnc::math::PiecewiseLinearFunction1d bound_line_1{{8.1, 11.9}, {60.0, 60.0}};
    tiny_pnc::math::PiecewiseLinearFunction1d bound_line_2{{3.1, 6.9}, {20.0, 20.0}};

    Motion1d motion;

    SUBCASE("SoftFences") {
        std::vector<SoftFence1d> soft_fences{
            {.id = 0,
             .actio = tiny_pnc::common::Actio::BLOCKING,
             .bound_line = bound_line_1,
             .linear_weight = 1000.0,
             .quadratic_weight = 10.0},
            {.id = 1,
             .actio = tiny_pnc::common::Actio::PUSHING,
             .bound_line = bound_line_2,
             .linear_weight = 1000.0,
             .quadratic_weight = 10.0}};

        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    SUBCASE("HardFences") {
        std::vector<HardFence1d> hard_fences{
            {.id = 0, .actio = tiny_pnc::common::Actio::BLOCKING, .bound_line = bound_line_1},
            {.id = 1, .actio = tiny_pnc::common::Actio::PUSHING, .bound_line = bound_line_2}};

        route_line_acc_model.setHardFences(hard_fences);
        motion = route_line_acc_model.solve().motion;

        CHECK_EQ(motion.s(bound_line_1.maxT()), doctest::Approx(60.0).epsilon(1E-2));
        CHECK_EQ(motion.s(bound_line_2.minT()), doctest::Approx(20.0).epsilon(1E-2));
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = tiny_pnc::math::linspace(0.0, 20.0, num_samples);
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
        {
            std::vector<axes_handle> ax;
            for (size_t i = 0; i < 5; ++i) {
                ax.emplace_back(subplot(2, 3, i));
            }
            ax[0]->grid(on);
            ax[0]->hold(on);
            ax[0]->plot(plot_ts, plot_ss, ".");
            ax[0]->plot(
                {bound_line_1.minT(), bound_line_1.maxT()},
                {bound_line_1.minY(), bound_line_1.maxY()}, "r-"
            );
            ax[0]->plot(
                {bound_line_1.minT(), bound_line_1.minT()}, {bound_line_1.minY(), 100.0}, "r-"
            );
            ax[0]->plot(
                {bound_line_1.maxT(), bound_line_1.maxT()}, {bound_line_1.minY(), 100.0}, "r-"
            );
            ax[0]->plot(
                {bound_line_2.minT(), bound_line_2.maxT()},
                {bound_line_2.minY(), bound_line_2.maxY()}, "r-"
            );
            ax[0]->plot(
                {bound_line_2.minT(), bound_line_2.minT()}, {0.0, bound_line_2.minY()}, "r-"
            );
            ax[0]->plot(
                {bound_line_2.maxT(), bound_line_2.maxT()}, {0.0, bound_line_2.maxY()}, "r-"
            );

            ax[1]->grid(on);
            ax[1]->plot(plot_ts, plot_vs, ".");

            ax[2]->grid(on);
            ax[2]->plot(plot_ts, plot_as, ".");

            ax[3]->grid(on);
            ax[3]->plot(plot_ts, plot_js, ".");

            ax[4]->grid(on);
            ax[4]->plot(plot_ts, plot_snaps, ".");
            show();

            for (std::size_t i = 0; i < 5; ++i) {
                ax[i]->clear();
            }
        }
    }
}

TEST_CASE("HeadwayScene") {
    constexpr std::size_t num_samples = 41;
    const std::vector<double> sample_ts = tiny_pnc::math::linspace(0.0, 20.0, num_samples);
    RouteLineAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 5.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(1.0);
    route_line_acc_model.setJerkCost(10.0);
    route_line_acc_model.setSnapCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    Motion1d motion;
    tiny_pnc::math::PiecewiseLinearFunction1d bound_line;

    SUBCASE("Yield") {
        bound_line = tiny_pnc::math::PiecewiseLinearFunction1d{{0.0, 20.0}, {-10.0, 40}};
        std::vector<SoftFence1d> soft_fences{
            {.id = 0,
             .actio = tiny_pnc::common::Actio::BLOCKING,
             .bound_line = bound_line,
             .linear_weight = 1000.0,
             .quadratic_weight = 10.0}};
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    SUBCASE("Overtake") {
        bound_line = tiny_pnc::math::PiecewiseLinearFunction1d{{0.0, 20.0}, {10.0, 210.0}};
        std::vector<SoftFence1d> soft_fences{
            {.id = 0,
             .actio = tiny_pnc::common::Actio::PUSHING,
             .bound_line = bound_line,
             .linear_weight = 1000.0,
             .quadratic_weight = 10.0}};
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = tiny_pnc::math::linspace(
            0.0 + tiny_pnc::math::kEpsilon, 20.0 - tiny_pnc::math::kEpsilon, num_samples
        );
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
        {
            std::vector<axes_handle> ax;
            for (std::size_t i = 0; i < 5; ++i) {
                ax.emplace_back(subplot(2, 3, i));
            }
            ax[0]->grid(on);
            ax[0]->hold(on);
            ax[0]->plot(
                {bound_line.minT(), bound_line.maxT()}, {bound_line.minY(), bound_line.maxY()}, "r-"
            );
            ax[0]->plot(plot_ts, plot_ss, ".");

            ax[1]->grid(on);
            ax[1]->plot(plot_ts, plot_vs, ".");

            ax[2]->grid(on);
            ax[2]->plot(plot_ts, plot_as, ".");

            ax[3]->grid(on);
            ax[3]->plot(plot_ts, plot_js, ".");

            ax[4]->grid(on);
            ax[4]->plot(plot_ts, plot_snaps, ".");
            show();

            for (std::size_t i = 0; i < 5; ++i) {
                ax[i]->clear();
            }
        }
    }
}

} // namespace kinetics
} // namespace tiny_pnc

int main(int argc, const char* argv[]) {
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
