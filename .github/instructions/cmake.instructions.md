---
applyTo: "CMakeLists.txt,**/*.cmake"
---

# CMake Coding Standards

## General Style

- Use CMake 3.31+ features
- Use **Ninja** as the generator
- Use **CMake Presets** for build configurations
- Indent with 2 or 4 spaces consistently (match surrounding code)

## Project CMake Functions

### Library Declaration

Use `boyle_cxx_library()` from `cmake/utils.cmake`:

```cmake
boyle_cxx_library(
    NAME library_name
    HDRS
        header1.hpp
        header2.hpp
    SRCS
        source.cpp
    DEPS
        dependency1
        dependency2
    PUBLIC
)
```

Options: `PUBLIC`, `TESTONLY`, `DISABLE_INSTALL`. Omitting `SRCS` creates an `INTERFACE` library.

### Module Declaration

Use `boyle_cxx_module()` for C++20 modules:

```cmake
boyle_cxx_module(
    NAME module_name
    IXXS
        module_interface.ixx
    DEPS
        dependency
    PUBLIC
)
```

### Test Declaration

Use `boyle_cxx_test()`:

```cmake
boyle_cxx_test(
    NAME component_test
    SRCS component_test.cpp
    DEPS component_library
)
```

Test targets automatically link `cxxopts`, `doctest`, and `Matplot++`.

## Dependency Management

- Dependencies managed via CPM in `cmake/third_party.cmake`
- CPM source cache: `${CMAKE_SOURCE_DIR}/third_party`

```cmake
CPMAddPackage(
    NAME package_name
    VERSION x.y.z
    GITHUB_REPOSITORY owner/repo
    OPTIONS
        "OPTION_NAME value"
)
```

## CMake Presets

Defined in `CMakePresets.json`:
- `linux-gcc-x64` - Linux GCC x86_64
- `linux-clang-x64` - Linux Clang x86_64
- `darwin-gcc-arm64` - macOS GCC ARM64
- `darwin-clang-arm64` - macOS Clang ARM64
- `windows-msvc-x64` - Windows MSVC x86_64

Developer-specific presets in `CMakeUserPresets.json` (symlinked from `developers/<name>/`).

## Best Practices

- Use `target_*` commands instead of global `add_*` commands
- Specify `PUBLIC`, `PRIVATE`, or `INTERFACE` for all `target_*` commands
- Prefer `boyle_cxx_library()` and `boyle_cxx_test()` over raw `add_library()`/`add_executable()`
- Use generator expressions for build/install interface differences
- Set `CMAKE_COMPILE_WARNING_AS_ERROR=ON` for all builds
- Export `compile_commands.json` for IDE/clang-tidy integration
