# C++ Coding Style

## Naming

| Element | Style | Example |
| ------- | ----- | ------- |
| Namespace | `snake_case` | `boyle::math` |
| Class / struct | `PascalCase` | `OsqpSolver` |
| Function / method | `camelCase` | `solveQP` |
| Public members | `snake_case` | `settings` |
| Private members | `m_` + `snake_case` | `m_workspace` |
| Constants | `kPascalCase` | `kEpsilon` |
| Template parameters | `PascalCase` | `Scalar`, `Index` |
| Macros | `BOYLE_UPPER_SNAKE` | `BOYLE_CHECK_PARAMS` |

## Trailing return types

**Required** for new work:

```cpp
auto norm(Vec3d const& v) noexcept -> double;
```

## Headers and includes

- `#pragma once`
- Order: standard library → third-party → `boyle/...`
- Never `using namespace` in headers

## Class design

- `final` when not meant for inheritance
- `explicit` unary constructors unless implicit conversion is intentional
- `override` on every virtual override
- `struct` for passive data; `class` for invariants

## GNU attributes

Apply on hot paths after profiling:

```cpp
[[using gnu: always_inline, leaf, hot]]
```

## Modern C++

- `constexpr` / `if constexpr` where possible
- Concepts for template constraints
- `std::span`, `std::string_view` for non-owning buffers/strings
- Prefer stack allocation and `alignas` for small fixed-size objects

## Namespaces

Close with a comment:

```cpp
} // namespace boyle::math
```

## Comments

- Doxygen file header (`@file`, `@brief`, `@copyright`) on new files
- Line comments only for non-obvious rationale, units, or tolerances

See `.github/instructions/cpp-coding-standards.instructions.md` for extended detail.
