/**
 * @file function1.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-09
 *
 * @copyright Copyright (c) 2023 Tiny-PnC Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <type_traits>
#include <vector>

#include "common/utils/macros.hpp"
#include "math/type_traits.hpp"

namespace tiny_pnc {
namespace math {

template <typename DerivedFunction1, typename T, typename U>
class [[nodiscard]] Function1 {
    static_assert(
        std::is_arithmetic_v<T> || isVecArithmeticV<T>,
        "The loaded type must has arithmetic operators."
    );
    static_assert(std::is_floating_point_v<U>, "The loaded type must be a floating-point type.");

  public:
    ENABLE_IMPLICIT_CONSTRUCTORS(Function1);
    virtual ~Function1() noexcept = default;

    [[using gnu: pure, always_inline]]
    T
    operator()(U t) const noexcept {
        return underlying()->operator()(t);
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t) const noexcept {
        return underlying()->derivative(t);
    }

    [[using gnu: pure, always_inline]]
    T derivative(U t, unsigned int order) const {
        return underlying()->derivative(t, order);
    }

    [[using gnu: pure, always_inline]]
    T integral(U lower_bound, U upper_bound) const noexcept {
        return underlying()->integral(lower_bound, upper_bound);
    }

    [[using gnu: pure, always_inline]]
    U minT() const noexcept {
        return underlying()->minT();
    }

    [[using gnu: pure, always_inline]]
    U maxT() const noexcept {
        return underlying()->maxT();
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T>
    operator()(const std::vector<U>& ts) const noexcept {
        const std::size_t size{ts.size()};
        std::vector<T> ys;
        ys.reserve(size);
        for (const U t : ts) {
            ys.emplace_back(operator()(t));
        }
        return ys;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> derivative(const std::vector<U>& ts) const noexcept {
        const std::size_t size{ts.size()};
        std::vector<T> dys;
        dys.reserve(size);
        for (const U t : ts) {
            dys.emplace_back(derivative(t));
        }
        return dys;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<T> derivative(const std::vector<U> ts, unsigned int order) const {
        const std::size_t size{ts.size()};
        std::vector<T> dnys;
        dnys.reserve(size);
        for (const U t : ts) {
            dnys.emplace_back(derivative(t, order));
        }
        return dnys;
    }

  private:
    [[using gnu: pure, always_inline]]
    DerivedFunction1* underlying() noexcept {
        return static_cast<DerivedFunction1*>(this);
    }

    [[using gnu: pure, always_inline]]
    const DerivedFunction1* underlying() const noexcept {
        return static_cast<const DerivedFunction1*>(this);
    }
};

} // namespace math
} // namespace tiny_pnc
