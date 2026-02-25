#pragma once

#if defined(_MSC_VER) && defined(__cplusplus)
#include <complex>
#pragma warning(disable : 4190)
#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#define lapack_complex_float_real(z) ((z).real())
#define lapack_complex_float_imag(z) ((z).imag())
#define lapack_complex_double_real(z) ((z).real())
#define lapack_complex_double_imag(z) ((z).imag())
#endif
