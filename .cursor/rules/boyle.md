# Boyle Library - Cursor AI Rules

> High-performance C++23 mathematical library for autonomous driving and robotics

## Project Context

This is the **Boyle Library**, providing fundamental numeric implementations for autonomous driving and robotics:

- **Piecewise polynomial functions/curves** (linear, cubic, quintic)
- **Sparse matrix formats** (COO, CSC, CSR, DOK, LIL)
- **Convex optimization** (BFGS, L-BFGS, Amoeba, OSQP)
- **Dense linear algebra** with optional BLAS/LAPACK acceleration

## Technical Requirements

| Requirement  | Specification                      |
| ------------ | ---------------------------------- |
| C++ Standard | C++23 with GNU extensions          |
| CMake        | 3.31+ with Ninja                   |
| Compilers    | GCC 14+ or Clang 19+               |
| Testing      | doctest 2.4.12                     |
| Formatting   | clang-format 19+ (Microsoft style) |
| Linting      | clang-tidy 19+                     |

## Code Style Rules

### Naming Conventions (MANDATORY)

```
Classes/Structs:    PascalCase + final     → Vector, PiecewiseCubicCurve, BfgsSolver
Template Params:    PascalCase             → Scalar, Size, NRows, NCols, Alloc
Type Aliases:       snake_case             → value_type, size_type, param_type
Functions:          camelCase              → eval(), euclidean(), selfConjugated()
Member Variables:   m_ prefix              → m_data, m_vec_of_s
Static Constants:   k prefix               → kSize, kDuplicateCriterion
Local Variables:    camelCase              → result, alpha, tempArcLengths
Files:              snake_case.hpp/.cpp    → piecewise_cubic_curve.hpp
Tests:              *_test.cpp             → vector_test.cpp
```

### Return Type Style (MANDATORY)

Always use **trailing return type** syntax:

```cpp
// ✅ CORRECT
constexpr auto size() noexcept -> size_type;
auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type>;

// ❌ WRONG
constexpr size_type size() noexcept;
detail::DenseNormTraitT<value_type> euclidean() const noexcept;
```

### Header File Structure

Every header MUST follow this structure:

```cpp
/**
 * @file filename.hpp
 * @author Author Name (email@example.com)
 * @brief Brief description of the file
 * @version 0.1
 * @date YYYY-MM-DD
 *
 * @copyright Copyright (c) YYYY Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

// 1. Standard library includes (alphabetical)
#include <algorithm>
#include <concepts>
#include <vector>

// 2. Third-party includes
#include "boost/serialization/access.hpp"

// 3. Project includes
#include "boyle/math/concepts.hpp"

namespace boyle::math {

// Implementation here

} // namespace boyle::math
```

### Class Declaration Pattern

```cpp
template <ScalarArithmetic Scalar, std::size_t Size>
class alignas(32) Vector final {
    friend class boost::serialization::access;

  public:
    // Type aliases
    using value_type = Scalar;
    using size_type = std::size_t;

    // Static constants
    static constexpr size_type kSize = Size;

    // Rule of five (explicit)
    constexpr Vector(const Vector& other) noexcept = default;
    constexpr auto operator=(const Vector& other) noexcept -> Vector& = default;
    constexpr Vector(Vector&& other) noexcept = default;
    constexpr auto operator=(Vector&& other) noexcept -> Vector& = default;
    constexpr ~Vector() noexcept = default;

    // Constructors
    [[using gnu: always_inline]]
    constexpr explicit Vector(size_type size) noexcept(!BOYLE_CHECK_PARAMS);

    // Member functions with GNU attributes
    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclidean() const noexcept -> value_type;

  private:
    // Serialization
    [[using gnu: always_inline]]
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void {
        archive & m_data;
    }

    std::array<value_type, Size> m_data;
};
```

### GNU Attributes (USE APPROPRIATELY)

```cpp
[[using gnu: pure]]              // No side effects, result depends only on args
[[using gnu: const]]             // Like pure but no pointer dereference
[[using gnu: always_inline]]     // Force inline
[[using gnu: hot]]               // Frequently called
[[using gnu: leaf]]              // Doesn't call back into caller's TU
[[using gnu: flatten]]           // Inline all calls within function

// Common combinations:
[[using gnu: pure, always_inline, hot]]           // Hot path getters
[[using gnu: const, always_inline, leaf]]         // Static constexpr getters
[[using gnu: always_inline]]                      // Constructors
```

### Conditional Parameter Checking

```cpp
[[using gnu: always_inline]]
constexpr explicit Vector(size_type size) noexcept(!BOYLE_CHECK_PARAMS) {
#if BOYLE_CHECK_PARAMS == 1
    if (size != kSize) [[unlikely]] {
        throw std::invalid_argument("Vector constructor: size mismatches!");
    }
#endif
}
```

### Move Semantics Pattern

Provide all four operator overloads for binary operations:

```cpp
// lvalue + lvalue
[[using gnu: pure, always_inline, hot]]
constexpr auto operator+(const Vector& obj) const& noexcept -> Vector;

// lvalue + rvalue (reuse rvalue)
[[using gnu: always_inline, hot]]
constexpr auto operator+(Vector&& obj) const& noexcept -> Vector&&;

// rvalue + lvalue (reuse this)
[[using gnu: always_inline, hot]]
constexpr auto operator+(const Vector& obj) && noexcept -> Vector&&;

// rvalue + rvalue (reuse this)
[[using gnu: always_inline, hot]]
constexpr auto operator+(Vector&& obj) && noexcept -> Vector&&;
```

### BLAS/LAPACK Integration Pattern

```cpp
[[using gnu: pure, always_inline, hot]]
constexpr auto euclidean() const noexcept -> value_type {
    value_type result(0.0);
#ifdef BOYLE_USE_BLAS_LAPACK
    if constexpr (std::is_same_v<value_type, float>) {
        result = cblas_snrm2(size(), data(), 1);
    } else if constexpr (std::is_same_v<value_type, double>) {
        result = cblas_dnrm2(size(), data(), 1);
    } else if constexpr (std::is_same_v<value_type, std::complex<float>>) {
        result = cblas_scnrm2(size(), data(), 1);
    } else if constexpr (std::is_same_v<value_type, std::complex<double>>) {
        result = cblas_dznrm2(size(), data(), 1);
    } else {
        // Fallback implementation
    }
#else
    // Standard library implementation
    result = std::sqrt(std::transform_reduce(...));
#endif
    return result;
}
```

## Namespace Structure

```
boyle                           # Root namespace
├── boyle::math                 # Mathematical primitives
│   ├── boyle::math::detail     # Implementation details (private)
│   └── boyle::math::pmr        # PMR allocator variants
├── boyle::cvxopm               # Convex optimization
├── boyle::kinetics             # Motion planning
└── boyle::common               # Common utilities
```

## Unit Test Pattern

```cpp
/**
 * @file component_test.cpp
 * @author Author Name (email@example.com)
 * @brief
 * @version 0.1
 * @date YYYY-MM-DD
 *
 * @copyright Copyright (c) YYYY Boyle Development Team
 *            All rights reserved.
 *
 */

#include "boyle/math/dense/vector.hpp"

#include <array>
#include <sstream>

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace {

constexpr std::size_t kNRows{16};
constexpr std::array<double, kNRows> kTestValues{...};

} // namespace

namespace boyle::math {

// clang-format off
TEST_CASE_TEMPLATE("VectorTest", T,
    Vector<float, 16>, Vector<double, 16>,
    VectorX<float>, VectorX<double>) {
    // clang-format on

    T vec(kNRows);
    // ... test implementation

    // Floating-point comparison with tolerance
    CHECK_EQ(actual, doctest::Approx(expected).epsilon(1E-8));

    // Serialization round-trip test
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << vec;

    T other_vec;
    std::istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> other_vec;

    CHECK_EQ(vec, other_vec);
}

} // namespace boyle::math
```

## CMake Patterns

### Library Declaration

```cmake
boyle_cxx_library(
    NAME vector
    HDRS
        vector.hpp
        vector_view.hpp
    DEPS
        Boost::headers
        Boost::serialization
    PUBLIC
)
```

### Test Declaration

```cmake
boyle_cxx_test(
    NAME vector_test
    SRCS vector_test.cpp
    DEPS vector
)
```

## Concepts (C++20)

Use project-defined concepts from `boyle/math/concepts.hpp`:

```cpp
// Arithmetic types
ScalarArithmetic    // std::is_arithmetic_v<T>
ComplexArithmetic   // std::complex<T>

// Linear algebra types
VecArithmetic       // Vector-like types with data(), size(), dot()
MatArithmetic       // Matrix-like types with nrows(), ncols()

// Combined
GeneralArithmetic   // ScalarArithmetic || VecArithmetic || MatArithmetic

// Other
Allocatory          // Allocator concept
Iterable            // Has begin()/end()
```

## Common Project Macros

From `boyle/common/utils/macros.hpp`:

```cpp
ENABLE_COPY(ClassName)
DISABLE_COPY(ClassName)
ENABLE_MOVE(ClassName)
DISABLE_MOVE(ClassName)
ENABLE_COPY_AND_MOVE(ClassName)
DISABLE_COPY_AND_MOVE(ClassName)
ENABLE_IMPLICIT_CONSTRUCTORS(ClassName)
DISABLE_IMPLICIT_CONSTRUCTORS(ClassName)
MAKE_SINGLETON(ClassName)
```

## Build Commands

```bash
# Configure
cmake --preset linux-gcc-x64

# Build
cmake --build --preset linux-gcc-x64

# Test
ctest --preset linux-gcc-x64

# Run specific test with visualization
./out/build/linux-gcc-x64/tests/boyle/math/curves/piecewise_cubic_curve_test --plot-graph
```

## Dependencies Quick Reference

| Library             | Header                           | Purpose              |
| ------------------- | -------------------------------- | -------------------- |
| Boost.Serialization | `boost/serialization/access.hpp` | Binary serialization |
| doctest             | `doctest/doctest.h`              | Unit testing         |
| spdlog              | `spdlog/spdlog.h`                | Logging              |
| cxxopts             | `cxxopts.hpp`                    | CLI parsing          |
| Matplot++           | `matplot/matplot.h`              | Test visualization   |
| OpenBLAS            | `cblas.h`, `lapacke.h`           | BLAS/LAPACK          |

## Code Quality Checklist

When generating code, ensure:

- [ ] Doxygen file header with author, date, copyright
- [ ] `#pragma once` guard
- [ ] Includes sorted: standard → third-party → project
- [ ] Trailing return type syntax
- [ ] `final` keyword on classes
- [ ] GNU attributes on member functions
- [ ] `constexpr` where applicable
- [ ] `noexcept` specification (with `!BOYLE_CHECK_PARAMS` if throwing)
- [ ] Move semantics overloads for operators
- [ ] `alignas(32)` for SIMD-friendly types
- [ ] `[[likely]]`/`[[unlikely]]` for branch hints
- [ ] Parameter validation with `#if BOYLE_CHECK_PARAMS == 1`
- [ ] Namespace closing comment `} // namespace boyle::math`
