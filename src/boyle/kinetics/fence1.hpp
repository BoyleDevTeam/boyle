/**
 * @file fence1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/vector.hpp"

#include "boyle/common/utils/macros.hpp"
#include "boyle/kinetics/dualism.hpp"

namespace boyle::kinetics {

template <std::floating_point T>
class [[nodiscard]] Fence1 {
    friend class boost::serialization::access;

  public:
    explicit Fence1(
        std::uint64_t c_id, ::boyle::kinetics::Actio c_actio, std::vector<T> c_bound_ts,
        std::vector<T> c_bound_ss
    ) noexcept
        : id{c_id}, actio{c_actio}, bound_ts{std::move(c_bound_ts)},
          bound_ss{std::move(c_bound_ss)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(Fence1);
    virtual ~Fence1() noexcept = 0;

    std::uint64_t id{std::numeric_limits<std::uint64_t>::quiet_NaN()};
    ::boyle::kinetics::Actio actio{::boyle::kinetics::Actio::BLOCKING};
    std::vector<T> bound_ts{};
    std::vector<T> bound_ss{};

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & id;
        archive & actio;
        archive & bound_ts;
        archive & bound_ss;
        return;
    }
};

template <std::floating_point T>
inline Fence1<T>::~Fence1() noexcept = default;

template <std::floating_point T>
class [[nodiscard]] HardFence1 final : public Fence1<T> {
    friend class boost::serialization::access;

  public:
    explicit HardFence1(
        std::uint64_t c_id, ::boyle::kinetics::Actio c_actio, std::vector<T> c_bound_ts,
        std::vector<T> c_bound_ss
    ) noexcept
        : Fence1<T>{c_id, c_actio, std::move(c_bound_ts), std::move(c_bound_ss)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(HardFence1);
    ~HardFence1() noexcept override = default;

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive& boost::serialization::base_object<Fence1<T>>(*this);
        return;
    }
};

template <std::floating_point T>
class [[nodiscard]] SoftFence1 final : public Fence1<T> {
    friend class boost::serialization::access;

  public:
    explicit SoftFence1(
        std::uint64_t c_id, ::boyle::kinetics::Actio c_actio, std::vector<T> c_bound_ts,
        std::vector<T> c_bound_ss, T c_linear_weight, T c_quadratic_weight
    ) noexcept
        : Fence1<T>{c_id, c_actio, std::move(c_bound_ts), std::move(c_bound_ss)},
          linear_weight{c_linear_weight}, quadratic_weight{c_quadratic_weight} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(SoftFence1);
    ~SoftFence1() noexcept override = default;

    T linear_weight{0.0};
    T quadratic_weight{0.0};

  private:
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive& boost::serialization::base_object<Fence1<T>>(*this);
        archive & linear_weight;
        archive & quadratic_weight;
        return;
    }
};

using HardFence1f = HardFence1<float>;
using HardFence1d = HardFence1<double>;

using SoftFence1f = SoftFence1<float>;
using SoftFence1d = SoftFence1<double>;

} // namespace boyle::kinetics
