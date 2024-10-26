/**
 * @file fft.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-10-22
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

#include <complex>
#include <concepts>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pocketfft_hdronly.h"

#include "boyle/math/concepts.hpp"

namespace boyle::math {

enum class FftNorm : std::uint8_t {
    BACKWARD = 0,
    ORTHO = 1,
    FORWARD = 2
};

template <ComplexArithmetic Complex>
[[using gnu: pure]] [[nodiscard]]
inline auto fft(std::span<const Complex> input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<Complex> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) [[unlikely]] {
        throw std::invalid_argument("boyle::math::fft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) [[unlikely]] {
        throw std::invalid_argument("boyle::math::fft(): invalid norm argument detected!");
    }
#endif

    using Scalar = typename Complex::value_type;
    Scalar fct = 1.0;
    switch (norm) {
    case FftNorm::BACKWARD:
        break;
    case FftNorm::ORTHO:
        fct /= std::sqrt(num_pnts);
        break;
    case FftNorm::FORWARD:
        fct /= num_pnts;
        break;
    default:
        std::unreachable();
    }

    std::vector<Complex> output(num_pnts);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(Complex)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(Complex)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2c(
        shape, strides_in, strides_out, axes, pocketfft::FORWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point Scalar>
[[using gnu: pure]] [[nodiscard]]
inline auto rfft(std::span<const Scalar> input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<Scalar>> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) [[unlikely]] {
        throw std::invalid_argument("boyle::math::rfft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) [[unlikely]] {
        throw std::invalid_argument("boyle::math::rfft(): invalid norm argument detected!");
    }
#endif

    using Complex = std::complex<Scalar>;
    Scalar fct = 1.0;
    switch (norm) {
    case FftNorm::BACKWARD:
        break;
    case FftNorm::ORTHO:
        fct /= std::sqrt(num_pnts / 2 + 1);
        break;
    case FftNorm::FORWARD:
        fct /= num_pnts / 2 + 1;
        break;
    default:
        std::unreachable();
    }

    std::vector<Complex> output(num_pnts / 2 + 1);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(Scalar)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(Complex)};
    const std::vector<std::size_t> axes{0};
    pocketfft::r2c(
        shape, strides_in, strides_out, axes, pocketfft::FORWARD, input.data(), output.data(), fct
    );

    return output;
}

template <ComplexArithmetic Complex>
[[using gnu: pure]] [[nodiscard]]
inline auto ifft(std::span<const Complex> input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<Complex> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) [[unlikely]] {
        throw std::invalid_argument("boyle::math::ifft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) [[unlikely]] {
        throw std::invalid_argument("boyle::math::ifft(): invalid norm argument detected!");
    }
#endif

    using Scalar = typename Complex::value_type;
    Scalar fct = 1.0;
    switch (norm) {
    case FftNorm::BACKWARD:
        fct /= num_pnts;
        break;
    case FftNorm::ORTHO:
        fct /= std::sqrt(num_pnts);
        break;
    case FftNorm::FORWARD:
        break;
    default:
        std::unreachable();
    }

    std::vector<Complex> output(num_pnts);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(Complex)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(Complex)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2c(
        shape, strides_in, strides_out, axes, pocketfft::BACKWARD, input.data(), output.data(), fct
    );

    return output;
}

template <ComplexArithmetic Complex>
[[using gnu: pure]] [[nodiscard]]
inline auto irfft(std::span<const Complex> input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<typename Complex::value_type> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) [[unlikely]] {
        throw std::invalid_argument("boyle::math::irfft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) [[unlikely]] {
        throw std::invalid_argument("boyle::math::irfft(): invalid norm argument detected!");
    }
#endif

    using Scalar = typename Complex::value_type;
    Scalar fct = 1.0;
    switch (norm) {
    case FftNorm::BACKWARD:
        fct /= num_pnts * 2 - 2;
        break;
    case FftNorm::ORTHO:
        fct /= std::sqrt(num_pnts * 2 - 2);
        break;
    case FftNorm::FORWARD:
        break;
    default:
        std::unreachable();
    }

    std::vector<Scalar> output(num_pnts * 2 - 2);

    const std::vector<std::size_t> shape{num_pnts * 2 - 2};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(Complex)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(Scalar)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2r(
        shape, strides_in, strides_out, axes, pocketfft::BACKWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point Scalar>
[[using gnu: const, always_inline]] [[nodiscard]]
inline auto fftfreq(std::size_t n, Scalar spacing = 1.0) noexcept -> std::vector<Scalar> {
    const Scalar val = 1.0 / (n * spacing);
    const int mid = (n - 1) / 2 + 1;
    const int n2 = n / 2;

    std::vector<Scalar> output(n);
    for (int i = 0; i < mid; ++i) {
        output[i] = i * val;
    }
    for (int i = mid; i < static_cast<int>(n); ++i) {
        output[i] = (i - mid - n2) * val;
    }

    return output;
}

template <std::floating_point Scalar>
[[using gnu: const, always_inline]] [[nodiscard]]
inline auto rfftfreq(std::size_t n, Scalar spacing = 1.0) noexcept -> std::vector<Scalar> {
    const Scalar val = 1.0 / (n * spacing);
    const std::size_t mid = n / 2 + 1;

    std::vector<Scalar> output(mid);
    for (std::size_t i = 0; i < mid; ++i) {
        output[i] = i * val;
    }

    return output;
}

} // namespace boyle::math
