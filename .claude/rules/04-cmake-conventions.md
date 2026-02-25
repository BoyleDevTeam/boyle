# CMake Conventions

## Helper functions

Defined in `cmake/utils.cmake`:

### `boyle_cxx_library`

- **NAME**: CMake target name (also `boyle::NAME` alias).
- **HDRS** / **SRCS**: Headers always listed; sources optional (interface library if only headers).
- **DEPS**: `PUBLIC` link libraries.
- Optional: `COPTS`, `DEFINES`, `LINKOPTS`, `PUBLIC`, `TESTONLY`, `DISABLE_INSTALL`.

### `boyle_cxx_module`

- For C++ module interfaces: **IXXS** file set, plus **SRCS**, **DEPS**, etc.

### `boyle_cxx_test`

- Builds a test executable when `BOYLE_BUILD_TESTING` is enabled, links doctest (and project-provided test deps), and registers `add_test`.

## CPM

- Third-party packages are fetched via **CPM** with pinned tags/commits.
- Keep dependency logic in dedicated `.cmake` fragments—do not scatter `CPMAddPackage` calls across leaf `CMakeLists.txt` unless the project already does so for that module.

## Presets

- Developers should use **`cmake --preset`** workflows from `CMakePresets.json`.
- Toolchain files live in `cmake/toolchains/`.

## Install

- Guard packaging rules with **`BOYLE_ENABLE_INSTALL`** to support CI and local dev without install.
