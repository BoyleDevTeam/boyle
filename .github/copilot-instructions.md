# Boyle — GitHub Copilot Instructions

This file orients coding agents and contributors to the **Boyle** codebase. Deeper, file-type-specific guidance lives under [`.github/instructions/`](instructions/).

## Project overview

Boyle is a **high-performance C++23** mathematical library (with **GNU extensions**) for autonomous driving and robotics: dense/sparse linear algebra, piecewise curves and functions, convex optimization, and motion-planning-oriented kinetics. The license is **BSD 3-Clause**. **Python** APIs and tooling are **planned** (uv, ruff, pyright, pytest).

## Technical stack

| Area | Choice |
| ---- | ------ |
| C++ | C++23, GCC 14+ / Clang 20+ / MSVC (Windows preset) |
| Build | CMake 3.31+ (Ninja), optional **xmake** |
| C++ deps | **CPM**, **vcpkg** (vendored tree) |
| Python deps | **uv** + `pyproject.toml` / `uv.lock` (when enabled) |
| Format | clang-format 20+ (Microsoft-based), ruff (Python) |
| Lint | clang-tidy 20+, ruff, pyright |
| Test | doctest 2.4.12 (C++), pytest (Python) |

## Module architecture

```
src/boyle/
├── common/    # Macros, logging, allocators, shared utilities
├── math/      # dense, sparse, curves, functions, mdfunctions
├── cvxopm/    # Convex optimization problems & solvers
└── kinetics/  # Motion models
```

Dependency order: **`common → math → cvxopm → kinetics`** (acyclic). Tests mirror this under `tests/boyle/`.

## Namespaces

- `boyle::common`, `boyle::math`, `boyle::math::pmr`, `boyle::math::detail`, `boyle::cvxopm`, `boyle::kinetics`
- Never `using namespace` in headers; close namespaces with `} // namespace boyle::...`

## Build commands

**CMake** (preferred):

```bash
cmake --preset linux-gcc-x64
cmake --build --preset linux-gcc-x64
ctest --preset linux-gcc-x64
ctest --preset linux-gcc-x64 -R <regex>
```

Presets: `linux-gcc-x64`, `linux-clang-x64`, `darwin-gcc-arm64`, `darwin-clang-arm64`, `windows-msvc-x64`. Build dir: `out/build/<preset>`.

**xmake**:

```bash
xmake build
```

**Python** (when present):

```bash
uv sync
uv run pytest
```

## Build options (CMake)

| CMake option | Purpose |
| ------------ | ------- |
| `BOYLE_CHECK_PARAMS` | Extra validation (`BOYLE_CHECK_PARAMS` macro) |
| `BOYLE_BUILD_TESTING` | Tests |
| `BOYLE_ENABLE_INSTALL` | Install rules |
| `BOYLE_USE_BLAS_LAPACK` | BLAS/LAPACK backend |
| `BOYLE_USE_BOOST_UNORDERED` | Boost unordered containers |
| `BOYLE_USE_SIMD` | SIMD kernels |

## Dependencies (indicative)

| Library | Role |
| ------- | ---- |
| Boost | Serialization, unordered |
| doctest | Unit tests |
| spdlog | Logging |
| OSQP / qdldl | QP |
| OpenBLAS | BLAS/LAPACK |
| cxxopts / Matplot++ | Test utilities |

Pins are in CMake/xmake dependency files—not duplicated here.

## Naming conventions

| Element | Style |
| ------- | ----- |
| Types | `PascalCase` |
| Functions | `camelCase` |
| Private fields | `m_snake_case` |
| Constants | `kPascalCase` |
| Macros | `BOYLE_SNAKE_CASE` |
| Type-alias suffix | `s` float, `d` double, `c` complex float, `z` complex double |

## Coding rules (summary)

1. **Trailing return type** on new functions: `auto foo() -> Bar;`
2. **`#pragma once`**; include order: std → third party → `boyle/...`
3. `final`, `explicit`, `override` per guidelines in [instructions](instructions/cpp-coding-standards.instructions.md)
4. **Exceptions** for errors; **`BOYLE_CHECK_PARAMS`** for optional defensive checks—not alternate status types
5. **Boost.Serialization** roundtrip tests for serializable types
6. **GNU attributes** on hot paths only (`[[using gnu: ...]]`)

## Error handling pattern

```cpp
constexpr auto coeff(size_type i) const noexcept(!BOYLE_CHECK_PARAMS) -> value_type {
#if BOYLE_CHECK_PARAMS == 1
  if (i >= kSize) [[unlikely]] {
    throw std::out_of_range("index out of range");
  }
#endif
  return m_data[i];
}
```

## Performance

- Stack + `alignas`; `std::span`; PMR where provided
- `constexpr` / concepts; avoid virtual/`std::function` in hot loops
- BLAS/LAPACK via configured backend; SIMD behind `BOYLE_USE_SIMD`

## Macros (`boyle/common/utils/macros.hpp`)

`ENABLE_COPY` / `DISABLE_COPY`, `ENABLE_MOVE` / `DISABLE_MOVE`, `ENABLE_COPY_AND_MOVE` / `DISABLE_COPY_AND_MOVE`, `MAKE_SINGLETON`, etc.—use these instead of ad-hoc `= delete` when matching project style.

## CMake targets

Helper functions in `cmake/utils.cmake`:

```cmake
boyle_cxx_library(NAME <target> HDRS "..." [SRCS "..."] DEPS ...)
boyle_cxx_test(NAME <target>_test SRCS "..._test.cpp" DEPS ...)
boyle_cxx_module(NAME <target> IXXS "..." DEPS ...)  # when modules are used
```

## Testing (doctest)

- `TEST_CASE_TEMPLATE` for multi-scalar coverage
- `#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN` **once** before `#include "doctest/doctest.h"`
- `doctest::Approx` with explicit epsilon for floats
- Header-under-test included first

## Security

- No hardcoded credentials; no secrets in tracked config

## TODO format

```
// TODO(github_id): description [#issue]
```

## Instruction files

| File | Scope |
| ---- | ----- |
| [workflow.instructions.md](instructions/workflow.instructions.md) | General workflow |
| [cpp-coding-standards.instructions.md](instructions/cpp-coding-standards.instructions.md) | `*.hpp`, `*.cpp` |
| [cpp-testing.instructions.md](instructions/cpp-testing.instructions.md) | `*_test.cpp` |
| [cmake.instructions.md](instructions/cmake.instructions.md) | CMake |
| [copyright.instructions.md](instructions/copyright.instructions.md) | Headers |
| [general.instructions.md](instructions/general.instructions.md) | All files |
| [security.instructions.md](instructions/security.instructions.md) | Secrets |
| [todo-format.instructions.md](instructions/todo-format.instructions.md) | TODOs |
| [python-coding-standards.instructions.md](instructions/python-coding-standards.instructions.md) | Python |
| [python-testing.instructions.md](instructions/python-testing.instructions.md) | Python tests |
