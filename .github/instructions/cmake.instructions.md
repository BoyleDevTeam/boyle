---
applyTo: "**/CMakeLists.txt,**/*.cmake"
---
# CMake instructions

## Version

Require CMake **3.31+** (see root `CMakeLists.txt`).

## Helper functions (`cmake/utils.cmake`)

- **`boyle_cxx_library(NAME ... HDRS ... [SRCS ...] DEPS ...)`** — Creates `boyle::NAME` alias; headers use target file sets.
- **`boyle_cxx_test(NAME ... SRCS ... DEPS ...)`** — Registers doctest executable + CTest entry when testing is enabled.
- **`boyle_cxx_module(NAME ... IXXS ... DEPS ...)`** — For C++ modules when used.

## CPM

- Centralize `CPMAddPackage` calls in shared includes; pin Git tags or commit hashes.
- Do not fetch dependencies ad hoc from leaf lists unless the repository already follows that pattern for the subdirectory.

## Presets

Developers run:

```bash
cmake --preset <name>
cmake --build --preset <name>
```

Prefer cache variables and toolchain files from `CMakePresets.json` over per-developer hardcoded paths.

## Style

- Use target-based APIs: `target_link_libraries`, `target_compile_definitions`, etc.
- Keep `CMakeLists.txt` readable: one logical target per block when possible.
