---
applyTo: "**/*_test.cpp"
---

# C++ Testing Guidelines

## Test Framework

- Use **doctest 2.4.12** framework
- Test files should end with `_test.cpp`
- Each test file must define `DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN`

## Test File Template

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
TEST_CASE_TEMPLATE("ComponentTest", T,
    Type1, Type2, Type3) {
    // clang-format on

    T instance(kNRows);
    CHECK_EQ(actual, doctest::Approx(expected).epsilon(1E-8));
}

} // namespace boyle::math
```

## Testing Patterns

### Template Tests

Use `TEST_CASE_TEMPLATE` for multi-type testing:
```cpp
TEST_CASE_TEMPLATE("VectorTest", T,
    Vector<float, 16>, Vector<double, 16>,
    VectorX<float>, VectorX<double>) {
    // ...
}
```

### Test Variations

Use `SUBCASE()` for test variations within a single test case.

### Floating-Point Comparison

```cpp
CHECK_EQ(actual, doctest::Approx(expected).epsilon(1E-8));
```

### Serialization Round-Trip

Always test Boost serialization round-trip:
```cpp
std::ostringstream oss;
boost::archive::binary_oarchive oa(oss);
oa << vec;

T other_vec;
std::istringstream iss(oss.str());
boost::archive::binary_iarchive ia(iss);
ia >> other_vec;

CHECK_EQ(vec, other_vec);
```

### Test Visualization

Support `--plot-graph` for visual tests with Matplot++.

## Assertions

- Use `CHECK_EQ`, `CHECK_NE`, `CHECK`, `CHECK_FALSE` for value checks
- Use `CHECK_THROWS`, `CHECK_THROWS_AS` for exception testing
- Use `REQUIRE_*` variants for fatal assertions that stop test execution

## CMakeLists.txt for Tests

```cmake
boyle_cxx_test(
    NAME component_test
    SRCS component_test.cpp
    DEPS component_library
)
```
