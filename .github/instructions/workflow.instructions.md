---
applyTo: "**/*"
---
# Workflow

## Development steps

1. Implement in the correct Boyle module (`common` → `math` → `cvxopm` → `kinetics`).
2. Add or update tests (`tests/boyle/` for C++; pytest when Python exists).
3. Format: **clang-format** on C++; **ruff** on Python.
4. Lint: **clang-tidy**; **pyright** + **ruff** for Python.
5. Build and run tests (CMake or xmake).

## File creation rules

- Do not add demos, tutorials, or unsolicited README/documentation files unless explicitly requested.
- Add C++ dependencies through **CPM** / **vcpkg** flows already in `cmake/`.
- Add Python dependencies with **`uv add`** and refresh **`uv lock`**.

## Build systems

- **Primary**: CMake + Ninja via presets (`CMakePresets.json`).
- **Alternative**: xmake (`xmake.lua`, `xmake/`).
- **Python**: `uv sync` for environments.

## Dependencies

- C++: pinned in CMake/CPM and vcpkg integration—follow existing patterns.
- Python: `pyproject.toml` + `uv.lock`.

## Performance hints

Prefer compile-time checks (`constexpr`, concepts), avoid heap in tight loops, use BLAS/LAPACK instead of naive linear algebra, and gate SIMD with `BOYLE_USE_SIMD`.
