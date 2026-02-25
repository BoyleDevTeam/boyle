---
applyTo: "**/*.hpp,**/*.cpp"
---

# C++ Coding Standards

## Naming Conventions

| Element          | Convention   | Examples                                              |
| ---------------- | ------------ | ----------------------------------------------------- |
| Classes/Structs  | PascalCase   | `Vector`, `PiecewiseCubicCurve`, `BfgsSolver`        |
| Template Params  | PascalCase   | `Scalar`, `Size`, `NRows`, `NCols`, `Alloc`          |
| Type Aliases     | snake_case   | `value_type`, `size_type`, `param_type`               |
| Concepts         | PascalCase   | `Arithmetic`, `ScalarArithmetic`, `VecArithmetic`     |
| Functions        | camelCase    | `eval()`, `euclidean()`, `selfConjugated()`           |
| Member Variables | `m_` prefix  | `m_data`, `m_vec_of_s`, `m_anchor_points`            |
| Static Constants | `k` prefix   | `kSize`, `kNRows`, `kDuplicateCriterion`              |
| Local Variables  | camelCase    | `result`, `alpha`, `ratio`                            |
| Enumerators      | `k` prefix   | `kSuccess`, `kFailure`                                |
| Headers          | snake_case   | `vector.hpp`, `piecewise_cubic_curve.hpp`             |
| Sources          | snake_case   | `osqp_solver.cpp`                                     |
| Tests            | `*_test.cpp` | `vector_test.cpp`, `matrix_test.cpp`                  |

## Trailing Return Types (MANDATORY)

Always use trailing return type syntax:
```cpp
constexpr auto size() noexcept -> size_type;
auto euclidean() const noexcept -> detail::DenseNormTraitT<value_type>;
```

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

## Headers and Includes

- Use `#pragma once` for header guards
- Include order: corresponding `.hpp` first, standard library (alphabetical), third-party, project headers
- Use `""` for third-party and project headers, `<>` for standard library only

## Class Design

- Mark all non-polymorphic classes with `final`
- Use `explicit` for single-argument constructors
- Use `struct` for passive data bundles; `class` for objects with invariants or behavior

### Class Declaration Pattern

```cpp
template <ScalarArithmetic Scalar, std::size_t Size>
class alignas(32) Vector final {
    friend class boost::serialization::access;

  public:
    using value_type = Scalar;
    static constexpr size_type kSize = Size;

    // Rule of five (explicit)
    constexpr Vector(const Vector& other) noexcept = default;
    constexpr auto operator=(const Vector& other) noexcept -> Vector& = default;

    [[using gnu: always_inline]]
    constexpr explicit Vector(size_type size) noexcept(!BOYLE_CHECK_PARAMS);

    [[using gnu: pure, always_inline, hot]]
    constexpr auto euclidean() const noexcept -> value_type;

  private:
    auto serialize(auto& archive, [[maybe_unused]] const unsigned int version) noexcept -> void;
    std::array<value_type, Size> m_data;
};
```

## GNU Attributes (USE APPROPRIATELY)

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
```

## Modern C++ Features

- Use C++20 concepts for template constraints (from `boyle/math/concepts.hpp`)
- Mark functions `constexpr` where applicable
- Use `noexcept` or `noexcept(!BOYLE_CHECK_PARAMS)` for functions that may throw parameter validation errors
- Use `alignas(32)` for SIMD-friendly types
- Use `[[likely]]` / `[[unlikely]]` for branch hints
- Prefer `std::ranges` and range-based algorithms
- Use `std::format` for string formatting

### Parameter Validation

```cpp
#if BOYLE_CHECK_PARAMS == 1
    if (size != kSize) [[unlikely]] {
        throw std::invalid_argument("Error message");
    }
#endif
```

### Move Semantics

Provide all four operator overloads for binary operations:
```cpp
auto operator+(const Vector& obj) const& noexcept -> Vector;
auto operator+(Vector&& obj) const& noexcept -> Vector&&;
auto operator+(const Vector& obj) && noexcept -> Vector&&;
auto operator+(Vector&& obj) && noexcept -> Vector&&;
```

## BLAS/LAPACK Integration

- Conditional compilation with `#ifdef BOYLE_USE_BLAS_LAPACK`
- Type-specific dispatch using `if constexpr`
- Always provide fallback standard library implementation

## Boost Serialization

- Friend declaration: `friend class boost::serialization::access;`
- Private serialize method: `auto serialize(auto& archive, ...) noexcept -> void;`

## Namespaces

- **Never use `using namespace xxx`**
- Use nested namespace syntax: `namespace boyle::math { ... }`
- Always close with comment: `} // namespace boyle::math`

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
