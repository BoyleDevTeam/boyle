/**
 * @file fft_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2024-10-22
 *
 * @copyright Copyright (c) 2024 Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/fft.hpp"

#include <complex>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle::math {

using namespace std::literals::complex_literals;

TEST_CASE_TEMPLATE("FftTest", T, float, double, long double) {
    std::vector<std::complex<T>> data{1.0, 2.0, 3.0, 4.0};
    std::vector<std::complex<T>> expected_fft_result{
        std::complex<T>(10.0), std::complex<T>(-2.0, 2.0), std::complex<T>(-2.0),
        std::complex<T>(-2.0, -2.0)
    };
    std::vector<std::complex<T>> result = fft(data);

    CHECK_EQ(result.size(), expected_fft_result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        CHECK_EQ(result[i].real(), doctest::Approx(expected_fft_result[i].real()).epsilon(1E-8));
        CHECK_EQ(result[i].imag(), doctest::Approx(expected_fft_result[i].imag()).epsilon(1E-8));
    }
}

TEST_CASE_TEMPLATE("RfftTest", T, float, double, long double) {
    std::vector<T> data{1.0, 2.0, 3.0, 4.0};
    std::vector<std::complex<T>> expected_fft_result{
        std::complex<T>(10.0), std::complex<T>(-2.0, 2.0), std::complex<T>(-2.0)
    };
    std::vector<std::complex<T>> result = rfft(data);

    CHECK_EQ(result.size(), expected_fft_result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        CHECK_EQ(result[i].real(), doctest::Approx(expected_fft_result[i].real()).epsilon(1E-8));
        CHECK_EQ(result[i].imag(), doctest::Approx(expected_fft_result[i].imag()).epsilon(1E-8));
    }
}

TEST_CASE_TEMPLATE("IfftTest", T, float, double, long double) {
    std::vector<std::complex<T>> fft_data{
        std::complex<T>(10.0), std::complex<T>(-2.0, 2.0), std::complex<T>(-2.0),
        std::complex<T>(-2.0, -2.0)
    };
    std::vector<std::complex<T>> expected_ifft_result{1.0, 2.0, 3.0, 4.0};
    std::vector<std::complex<T>> result = ifft(fft_data);

    CHECK_EQ(result.size(), expected_ifft_result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        CHECK_EQ(result[i].real(), doctest::Approx(expected_ifft_result[i].real()).epsilon(1E-8));
        CHECK_EQ(result[i].imag(), doctest::Approx(expected_ifft_result[i].imag()).epsilon(1E-8));
    }
}

TEST_CASE_TEMPLATE("IrfftTest", T, float, double, long double) {
    std::vector<std::complex<T>> rfft_data{
        std::complex<T>(10.0), std::complex<T>(-2.0, 2.0), std::complex<T>(-2.0)
    };
    std::vector<T> expected_ifft_result{1.0, 2.0, 3.0, 4.0};
    std::vector<T> result = irfft(rfft_data);

    CHECK_EQ(result.size(), expected_ifft_result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        CHECK_EQ(result[i], doctest::Approx(expected_ifft_result[i]).epsilon(1E-8));
    }
}

TEST_CASE_TEMPLATE("FftfreqTest", T, float, double, long double) {
    std::vector<T> expected_freq{0.0, 1.25, 2.5, 3.75, -5.0, -3.75, -2.5, -1.25};
    std::vector<T> freq = fftfreq<T>(8, 0.1);

    CHECK_EQ(freq.size(), expected_freq.size());

    for (size_t i = 0; i < freq.size(); ++i) {
        CHECK_EQ(freq[i], doctest::Approx(expected_freq[i]).epsilon(1E-8));
    }
}

TEST_CASE_TEMPLATE("RfftfreqTest", T, float, double, long double) {
    std::vector<T> expected_freq{0.0, 1.25, 2.5, 3.75, 5.0};
    std::vector<T> freq = rfftfreq<T>(8, 0.1);

    CHECK_EQ(freq.size(), expected_freq.size());

    for (size_t i = 0; i < freq.size(); ++i) {
        CHECK_EQ(freq[i], doctest::Approx(expected_freq[i]).epsilon(1E-8));
    }
}

} // namespace boyle::math
