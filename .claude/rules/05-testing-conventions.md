# Testing Conventions (C++)

## Framework

- **doctest** 2.4.12

## Main macro

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
```

Place the `#define` **before** the doctest include. **Exactly one** TU per test binary should provide this.

## Templates

Use **`TEST_CASE_TEMPLATE`** to cover `float`, `double`, and complex scalar types where algorithms are generic.

## Includes

1. Header under test (when applicable)
2. Project headers
3. doctest

## Assertions

- **`doctest::Approx`** for floating comparisons with explicit tolerances:

```cpp
REQUIRE(x == doctest::Approx(y).epsilon(1e-12));
```

## Serialization

For Boost-serializable types, verify **serialize → deserialize → equal** to the original object (or approximate equality for floats).

## Naming

- Test executables and `TEST_CASE` strings should read like specifications (`"LU decomposition recovers identity on well-conditioned matrix"`).

## Organization

- Mirror `src/boyle/` layout under `tests/boyle/`.
- Avoid global state and order dependencies between tests.
