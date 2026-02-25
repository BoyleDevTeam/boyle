---
applyTo: "**/*.cpp,**/*.hpp"
---
# C++ coding standards

## Naming

- **Namespaces**: `snake_case` (`boyle::math`).
- **Classes / structs**: `PascalCase`.
- **Functions / methods**: `camelCase`.
- **Private members**: `m_` prefix + `snake_case`.
- **Constants**: `kPascalCase`.
- **Macros**: `BOYLE_UPPER_SNAKE`.

## Trailing return types

Use for new and updated functions:

```cpp
auto frobnicate(Vec3d x) noexcept -> double;
```

## Headers

- `#pragma once`.
- Include order: C++ standard library, third party, then `boyle/...` headers.
- No `using namespace` in headers.

## Class design

- `final` for non-polymorphic types; `explicit` constructors unless implicit conversion is intended; `override` on overrides; virtual destructor when base is polymorphic.
- Use project macros (`ENABLE_COPY`, `DISABLE_MOVE`, …) from `boyle/common/utils/macros.hpp` when appropriate.

## GNU attributes

Apply on proven hot paths:

```cpp
[[using gnu: always_inline, leaf, hot]]
```

## Modern C++

- Concepts for template constraints; `if constexpr` where it simplifies generics.
- `std::span` / `std::string_view` for non-owning data.
- Exceptions for unrecoverable misuse; `BOYLE_CHECK_PARAMS` for optional validation.

## Serialization

Use Boost.Serialization patterns already present in the module; add roundtrip tests.
