# Boyle Library - GitHub Copilot Instructions

## Project Overview

Boyle Library is a high-performance C++23 mathematical library for autonomous driving and robotics. It provides piecewise polynomial functions/curves, sparse matrices, and convex optimization modules for trajectory generation.

## Technical Stack

- **Language**: C++23 with GNU extensions
- **Build System**: CMake 3.31+ with Ninja generator
- **Package Manager**: CPM (CMake Package Manager)
- **Compilers**: GCC 14+ or Clang 19+
- **Testing**: doctest 2.4.12
- **Static Analysis**: clang-tidy 19+, clang-format 19+

## Code Architecture

### Module Structure

```
src/boyle/
├── common/       # Utilities: FSM, allocators, logging, macros
├── math/         # Core math: vectors, matrices, curves, functions, FFT
│   ├── dense/    # Dense linear algebra (Vector, Matrix, LU, QR)
│   ├── sparse/   # Sparse matrices (COO, CSC, CSR, DOK, LIL)
│   ├── curves/   # Piecewise polynomial curves (linear, cubic, quintic)
│   ├── functions/# Piecewise polynomial functions
│   └── mdfunctions/ # Multi-dimensional functions
├── cvxopm/       # Convex optimization (BFGS, L-BFGS, OSQP)
└── kinetics/     # Motion planning models
```

### Namespace Organization

- `boyle` - Root namespace
- `boyle::math` - Mathematical primitives
- `boyle::math::pmr` - PMR allocator variants
- `boyle::math::detail` - Implementation details
- `boyle::cvxopm` - Optimization solvers
- `boyle::kinetics` - Motion planning
- `boyle::common` - Common utilities

## Naming Conventions

### Classes and Types

- **Classes**: PascalCase with `final` keyword - `Vector`, `Matrix`, `PiecewiseCubicCurve`, `BfgsSolver`
- **Template Parameters**: PascalCase - `Scalar`, `Size`, `NRows`, `NCols`, `Order`, `Alloc`
- **Type Aliases**: snake_case - `value_type`, `size_type`, `allocator_type`, `param_type`
- **Concepts**: PascalCase - `ScalarArithmetic`, `VecArithmetic`, `MatArithmetic`, `Allocatory`

### Functions and Variables

- **Functions**: camelCase - `eval()`, `euclidean()`, `normalized()`, `setIdentity()`, `selfConjugated()`
- **Member Variables**: `m_` prefix - `m_data`, `m_vec_of_s`, `m_anchor_points`
- **Static Constants**: `k` prefix - `kSize`, `kNRows`, `kDuplicateCriterion`, `kFactors`
- **Local Variables**: camelCase - `result`, `alpha`, `ratio`, `temp_arc_lengths`

### Files

- **Headers**: snake_case with `.hpp` - `vector.hpp`, `piecewise_cubic_curve.hpp`
- **Sources**: snake_case with `.cpp` - `osqp_solver.cpp`
- **Tests**: `*_test.cpp` - `vector_test.cpp`, `matrix_test.cpp`

## Header File Template

```cpp
/**
 * @file filename.hpp
 * @author Author Name (email@example.com)
 * @brief Brief description
 * @version 0.1
 * @date YYYY-MM-DD
 *
 * @copyright Copyright (c) YYYY Boyle Development Team
 *            All rights reserved.
 *
 */

#pragma once

// Standard library includes (alphabetical)
#include <algorithm>
#include <concepts>
#include <vector>

// Third-party includes
#include "boost/serialization/access.hpp"

// Project includes
#include "boyle/math/concepts.hpp"

namespace boyle::math {

// Implementation

} // namespace boyle::math
```

## Modern C++ Patterns

### Required Features

1. **C++20 Concepts**: Use concepts for template constraints

    ```cpp
    template <ScalarArithmetic Scalar, std::size_t Size>
    class Vector final { ... };
    ```

2. **GNU Attributes**: Use `[[using gnu: ...]]` for optimization hints

    ```cpp
    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclidean() const noexcept -> ... { }
    ```

3. **Trailing Return Types**: Always use trailing return type syntax

    ```cpp
    constexpr auto operator[](size_type i) noexcept -> reference;
    ```

4. **Constexpr Functions**: Mark functions `constexpr` where applicable

5. **noexcept Specification**: Use `noexcept` or `noexcept(!BOYLE_CHECK_PARAMS)`

6. **Move Semantics**: Provide rvalue overloads for operators

    ```cpp
    auto operator+(const Vector& obj) const& noexcept -> Vector;
    auto operator+(Vector&& obj) const& noexcept -> Vector&&;
    ```

7. **Alignment**: Use `alignas(32)` for SIMD-friendly types

8. **Branch Hints**: Use `[[likely]]` and `[[unlikely]]`

### BLAS/LAPACK Integration

- Conditional compilation with `#ifdef BOYLE_USE_BLAS_LAPACK`
- Type-specific dispatch using `if constexpr`
- Fallback to standard library algorithms

### Boost Serialization

- Friend declaration: `friend class boost::serialization::access;`
- Private serialize method:
    ```cpp
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void;
    ```

## Unit Testing

### Test File Template

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

#include "boyle/path/to/component.hpp"

#include <array>
// other includes...

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace {

// Test constants
constexpr std::size_t kNRows{16};
constexpr std::array<double, kNRows> kTestValues{...};

} // namespace

namespace boyle::math {

TEST_CASE_TEMPLATE("ComponentTest", T, Type1, Type2, ...) {
    // Test implementation
    CHECK_EQ(actual, doctest::Approx(expected).epsilon(1E-8));
}

} // namespace boyle::math
```

### Testing Patterns

- Use `TEST_CASE_TEMPLATE` for multi-type testing
- Use `SUBCASE()` for test variations
- Use `doctest::Approx().epsilon()` for floating-point comparisons
- Test Boost serialization round-trip
- Support `--plot-graph` for visual tests with Matplot++

## CMake Integration

### Library Declaration

```cmake
boyle_cxx_library(
    NAME library_name
    HDRS header1.hpp header2.hpp
    SRCS source.cpp  # Optional for interface libraries
    DEPS dependency1 dependency2
    PUBLIC  # or TESTONLY
)
```

### Test Declaration

```cmake
boyle_cxx_test(
    NAME test_name
    SRCS test_file.cpp
    DEPS library_dependency
)
```

## Code Quality Requirements

1. **Warnings as Errors**: `CMAKE_COMPILE_WARNING_AS_ERROR=ON`
2. **clang-format**: Microsoft style base, 100 column limit, 4-space indent
3. **clang-tidy**: Comprehensive checks (bugprone, modernize, performance, readability, cppcoreguidelines)

## Dependencies

| Library   | Version | Purpose                |
| --------- | ------- | ---------------------- |
| Boost     | 1.90.0  | Serialization, headers |
| doctest   | 2.4.12  | Unit testing           |
| spdlog    | 1.16.0  | Logging                |
| cxxopts   | 3.3.1   | CLI parsing            |
| Matplot++ | 1.2.2   | Test visualization     |
| OSQP      | 1.0.0   | QP solving             |
| OpenBLAS  | -       | BLAS/LAPACK (optional) |

## Common Macros

```cpp
// From boyle/common/utils/macros.hpp
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

## Parameter Validation

Use conditional exception throwing:

```cpp
#if BOYLE_CHECK_PARAMS == 1
    if (size != kSize) [[unlikely]] {
        throw std::invalid_argument("Error message");
    }
#endif
```

## Hints for Code Generation

1. Always include Doxygen file headers
2. Prefer `std::ranges` and range-based algorithms
3. Use `std::format` for string formatting (C++20)
4. Prefer `std::span` over raw pointer + size
5. Use structured bindings where appropriate
6. Mark single-argument constructors `explicit`
7. Use `[[nodiscard]]` for functions returning important values
8. Avoid raw `new`/`delete` - use smart pointers or containers
