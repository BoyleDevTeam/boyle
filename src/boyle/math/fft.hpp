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

#include <algorithm>
#include <complex>
#include <concepts>
#include <cstdint>
#include <stdexcept>
#include <vector>

#define POCKETFFT_NO_MULTITHREADING
#include "pocketfft_hdronly.h"

namespace boyle::math {

enum class FftNorm : std::uint8_t {
    BACKWARD = 0,
    ORTHO = 1,
    FORWARD = 2
};

template <std::floating_point T>
[[using gnu: pure]] [[nodiscard]]
inline auto fft(const std::vector<std::complex<T>>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<T>> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) {
        throw std::invalid_argument("boyle::math::fft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) {
        throw std::invalid_argument("boyle::math::fft(): invalid norm argument detected!");
    }
#endif

    T fct{static_cast<T>(1.0)};
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
        break;
    }

    std::vector<std::complex<T>> output(num_pnts);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(std::complex<T>)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(std::complex<T>)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2c(
        shape, strides_in, strides_out, axes, pocketfft::FORWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]] [[nodiscard]]
inline auto fft(const std::vector<T>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<T>> {
    std::vector<std::complex<T>> complex_input(input.size());
    std::transform(
        input.cbegin(), input.cend(), complex_input.begin(),
        [](const T& x) noexcept -> std::complex<T> { return std::complex<T>(x, 0.0); }
    );
    return fft(complex_input, norm);
}

template <std::floating_point T>
[[using gnu: pure]] [[nodiscard]]
inline auto rfft(const std::vector<T>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<T>> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) {
        throw std::invalid_argument("boyle::math::rfft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) {
        throw std::invalid_argument("boyle::math::rfft(): invalid norm argument detected!");
    }
#endif

    T fct{static_cast<T>(1.0)};
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
        break;
    }

    std::vector<std::complex<T>> output(num_pnts / 2 + 1);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(T)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(std::complex<T>)};
    const std::vector<std::size_t> axes{0};
    pocketfft::r2c(
        shape, strides_in, strides_out, axes, pocketfft::FORWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point T>
[[using gnu: pure]] [[nodiscard]]
inline auto ifft(const std::vector<std::complex<T>>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<T>> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) {
        throw std::invalid_argument("boyle::math::ifft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) {
        throw std::invalid_argument("boyle::math::ifft(): invalid norm argument detected!");
    }
#endif

    T fct{static_cast<T>(1.0)};
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
        break;
    }

    std::vector<std::complex<T>> output(num_pnts);

    const std::vector<std::size_t> shape{num_pnts};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(std::complex<T>)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(std::complex<T>)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2c(
        shape, strides_in, strides_out, axes, pocketfft::BACKWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point T>
[[using gnu: pure, always_inline]] [[nodiscard]]
inline auto ifft(const std::vector<T>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<std::complex<T>> {
    std::vector<std::complex<T>> complex_input(input.size());
    std::transform(
        input.cbegin(), input.cend(), complex_input.begin(),
        [](const T& x) noexcept -> std::complex<T> { return std::complex<T>(x, 0.0); }
    );
    return ifft(complex_input, norm);
}

template <std::floating_point T>
[[using gnu: pure]] [[nodiscard]]
inline auto irfft(const std::vector<std::complex<T>>& input, FftNorm norm = FftNorm::BACKWARD)
    -> std::vector<T> {
    const std::size_t num_pnts{input.size()};
#if BOYLE_CHECK_PARAMS == 1
    if (num_pnts < 2) {
        throw std::invalid_argument("boyle::math::irfft(): input size must be greater than 1!");
    }
    if (norm > FftNorm::FORWARD) {
        throw std::invalid_argument("boyle::math::irfft(): invalid norm argument detected!");
    }
#endif

    T fct{static_cast<T>(1.0)};
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
        break;
    }

    std::vector<T> output(num_pnts * 2 - 2);

    const std::vector<std::size_t> shape{num_pnts * 2 - 2};
    const std::vector<std::ptrdiff_t> strides_in{sizeof(std::complex<T>)};
    const std::vector<std::ptrdiff_t> strides_out{sizeof(T)};
    const std::vector<std::size_t> axes{0};
    pocketfft::c2r(
        shape, strides_in, strides_out, axes, pocketfft::BACKWARD, input.data(), output.data(), fct
    );

    return output;
}

template <std::floating_point T>
[[using gnu: const, always_inline]] [[nodiscard]]
inline auto fftfreq(const std::size_t n, const T spacing = 1.0) noexcept -> std::vector<T> {
    const T val = static_cast<T>(1.0) / (n * spacing);
    const int mid = (n - 1) / 2 + 1;
    const int n2 = n / 2;

    std::vector<T> output(n);
    for (int i = 0; i < mid; ++i) {
        output[i] = i * val;
    }
    for (int i = mid; i < static_cast<int>(n); ++i) {
        output[i] = (i - mid - n2) * val;
    }

    return output;
}

template <std::floating_point T>
[[using gnu: const, always_inline]] [[nodiscard]]
inline auto rfftfreq(const std::size_t n, const T spacing = 1.0) noexcept -> std::vector<T> {
    const T val = static_cast<T>(1.0) / (n * spacing);
    const std::size_t mid = n / 2 + 1;

    std::vector<T> output(mid);
    for (std::size_t i = 0; i < mid; ++i) {
        output[i] = i * val;
    }

    return output;
}

} // namespace boyle::math
