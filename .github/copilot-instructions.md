# Boyle Library - GitHub Copilot Instructions

## Project Overview

Boyle Library is a high-performance C++23 mathematical library for autonomous driving and robotics. It provides piecewise polynomial functions/curves, sparse matrices, and convex optimization modules for trajectory generation.

## Technical Stack

- **Language**: C++23 with GNU extensions
- **Build System**: CMake 3.31+ with Ninja generator
- **Package Manager**: CPM (CMake Package Manager)
- **Compilers**: GCC 14+ or Clang 20+
- **Testing**: doctest 2.4.12
- **Static Analysis**: clang-tidy 20+, clang-format 20+ (Microsoft style)

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
├── cvxopm/       # Convex optimization (BFGS, L-BFGS, Amoeba, OSQP)
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

## Build Commands

```bash
cmake --preset linux-gcc-x64
cmake --build --preset linux-gcc-x64
ctest --preset linux-gcc-x64
```

## Build Configuration Options

| Option                       | Default    | Description                       |
| ---------------------------- | ---------- | --------------------------------- |
| `BOYLE_CHECK_PARAMS`         | `OFF`      | Enable parameter checking         |
| `BOYLE_BUILD_TESTING`        | `ON`       | Enable testing                    |
| `BOYLE_ENABLE_INSTALL`       | `ON`       | Enable install targets            |
| `BOYLE_USE_BLAS_LAPACK`      | `OpenBLAS` | BLAS/LAPACK backend               |
| `BOYLE_USE_BOOST_UNORDERED`  | `ON`       | Enable Boost.unordered containers |
| `BOYLE_USE_SIMD`             | `OFF`      | SIMD implementation (AVX512/OFF)  |
| `BUILD_SHARED_LIBS`          | `OFF`      | Build shared libraries            |
| `CMAKE_UNITY_BUILD`          | `OFF`      | Enable unity build                |

## Dependency Management

Dependencies managed via CPM in `cmake/third_party.cmake`. Use `CPMAddPackage()` to add new dependencies.

| Library         | Version | Purpose                  |
| --------------- | ------- | ------------------------ |
| Boost           | 1.90.0  | Serialization, headers   |
| doctest         | 2.4.12  | Unit testing             |
| spdlog          | 1.17.0  | Logging                  |
| cxxopts         | 3.3.1   | CLI parsing              |
| Microsoft Proxy | 4.0.2   | Type-erased polymorphism |
| Matplot++       | 1.2.2   | Test visualization       |
| OSQP            | 1.0.0   | QP solving               |
| OpenBLAS        | -       | BLAS/LAPACK (optional)   |

## Security Requirements

- **NEVER** hardcode API keys, tokens, or passwords in code
- **NEVER** commit `.env` files with real credentials
- Use environment variables for secrets

## TODO Format

```cpp
// TODO(github_id): description [optional issue reference]
```

## Code Quality Requirements

1. **Warnings as Errors**: `CMAKE_COMPILE_WARNING_AS_ERROR=ON`
2. **clang-format**: Microsoft style base, 100 column limit, 4-space indent, no tabs, pointer alignment left, attach braces
3. **clang-tidy**: Comprehensive checks (bugprone, cert, hicpp, misc, modernize, performance, readability, cppcoreguidelines, concurrency, portability, boost)
4. **Comments**: English only, don't comment obvious code
5. **File headers**: Doxygen format with `@file`, `@author`, `@brief`, `@version`, `@date`, `@copyright`

## Hints for Code Generation

1. Always include Doxygen file headers
2. Prefer `std::ranges` and range-based algorithms
3. Use `std::format` for string formatting (C++20)
4. Prefer `std::span` over raw pointer + size
5. Use structured bindings where appropriate
6. Mark single-argument constructors `explicit`
7. Use `[[nodiscard]]` for functions returning important values
8. Avoid raw `new`/`delete` - use smart pointers or containers
9. **NEVER** create demo programs, example files, or documentation unless explicitly requested
10. Run `clang-format -i {files}` after any code changes
