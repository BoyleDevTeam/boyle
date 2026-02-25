# AGENTS.md — Boyle Project Guide for AI Agents

## Project Overview

**Boyle** is a high-performance **C++23** mathematical library for autonomous driving and robotics. It provides dense and sparse linear algebra, piecewise curves and functions, convex optimization solvers, and motion-planning oriented kinetics. The project is licensed under the **BSD 3-Clause** license. **Python bindings and tooling are planned**; treat Python sections of this guide as forward-looking where noted.

The codebase uses **GNU extensions** where beneficial and targets multiple host platforms via CMake presets and an optional **xmake** workflow.

## Build System

- **Primary**: **CMake 3.31+** with **Ninja** (see `CMakePresets.json`).
- **Alternatives**: **xmake** (see root `xmake.lua` and `xmake/`).
- **C++ dependencies**: **CPM** (CMake Package Manager) and **vcpkg** (vendored under `vcpkg/` as applicable).
- **Python (future)**: **uv** with `pyproject.toml` / `uv.lock`.

### Key Commands

| Command                              | Purpose                                               |
| ------------------------------------ | ----------------------------------------------------- |
| `cmake --preset <preset>`            | Configure (e.g. `linux-gcc-x64`, `linux-clang-x64`)   |
| `cmake --build --preset <preset>`    | Build all targets                                     |
| `ctest --preset <preset>`            | Run all tests for that preset                         |
| `ctest --preset <preset> -R <regex>` | Run tests matching a name                             |
| `clang-format -i <files>`            | Format C/C++ (Microsoft-based style, see repo config) |
| `clang-tidy <args>`                  | C++ static analysis                                   |
| `xmake` / `xmake build`              | Configure and build via xmake                         |
| `uv sync` / `uv run pytest`          | Python env and tests (when Python is enabled)         |
| `ruff check` / `ruff format`         | Python lint and format (when Python is present)       |

### CMake Presets

| Preset               | Host                  |
| -------------------- | --------------------- |
| `linux-gcc-x64`      | Linux, GCC, x86_64    |
| `linux-clang-x64`    | Linux, Clang, x86_64  |
| `darwin-gcc-arm64`   | macOS, GCC, arm64     |
| `darwin-clang-arm64` | macOS, Clang, arm64   |
| `windows-msvc-x64`   | Windows, MSVC, x86_64 |

Configure and build output directories follow `CMakePresets.json` (typically `out/build/<preset>` and `out/install/<preset>`).

## Languages & Toolchains

| Language   | Version / Toolchain                                    | Primary Location                                     |
| ---------- | ------------------------------------------------------ | ---------------------------------------------------- |
| **C++**    | C++23, GCC 14+ or Clang 20+, MSVC where preset applies | `src/boyle/`                                         |
| **CMake**  | 3.31+                                                  | `CMakeLists.txt`, `cmake/`                           |
| **xmake**  | xmake 2.9+ (see `set_xmakever`)                        | `xmake.lua`, `xmake/`                                |
| **Python** | Future: 3.12+ via **uv**                               | (to be defined under `src/` or `python/` when added) |

## Code Formatting

- **C++**: **clang-format 20+**, **Microsoft** base style (see `.clang-format`). Format changed files before completing work.
- **Python (future)**: **ruff** for lint + format; **pyright** for type checking.
- **CMake**: follow project conventions and keep lists readable; no separate clang-format requirement unless configured.

## Repository Structure

```
boyle/
├── src/boyle/           # Library sources (headers + occasional .cpp)
│   ├── common/          # Utilities, allocators, logging, macros
│   ├── math/            # dense, sparse, curves, functions, mdfunctions
│   ├── cvxopm/          # Convex optimization (problems, solvers)
│   └── kinetics/        # Motion planning models
├── tests/boyle/         # doctest-based tests (mirror of src layout)
├── cmake/               # Toolchains, utils.cmake (boyle_cxx_* helpers), third_party hooks
├── xmake/               # xmake rules and third-party recipes
├── out/                 # Build trees (CMake/xmake; typically gitignored)
├── CMakePresets.json
├── CMakeLists.txt
├── xmake.lua
├── pyproject.toml       # Python metadata (future)
└── vcpkg/               # Vendored vcpkg (large)
```

## Module Architecture

Dependency direction (no cycles):

```
common → math → cvxopm → kinetics
```

- **common**: FSM helpers, allocators, logging (e.g. spdlog), shared macros.
- **math**: Vectors, matrices, sparse formats, curves, piecewise functions, multi-dimensional functions.
- **cvxopm**: QP and related convex optimization (OSQP, L-BFGS, etc., per `cmake/` / CPM).
- **kinetics**: Models for motion planning.

Library code is mostly header-only with fine-grained CMake (and xmake) targets per component.

## Namespace Structure

| Namespace             | Role                                |
| --------------------- | ----------------------------------- |
| `boyle`               | Root                                |
| `boyle::common`       | Shared utilities                    |
| `boyle::math`         | Core math types and algorithms      |
| `boyle::math::pmr`    | Polymorphic allocator variants      |
| `boyle::math::detail` | Internal implementation details     |
| `boyle::cvxopm`       | Optimization interfaces and solvers |
| `boyle::kinetics`     | Kinematic / planning models         |

**Never** use `using namespace` at namespace scope in headers. Close namespaces with a trailing comment (e.g. `} // namespace boyle::math`).

## Development Workflow

1. **Implement** functionality in the appropriate module under `src/boyle/`.
2. **Add or update tests** under `tests/boyle/` (doctest).
3. **Build** with `cmake --preset ...` and `cmake --build --preset ...`, or xmake.
4. **Run tests**: `ctest --preset ...` (or xmake test targets when using xmake).
5. **Format**: `clang-format -i` on touched C++ files; **ruff** on Python when applicable.
6. **Lint**: `clang-tidy` on C++; **ruff** + **pyright** for Python when present.
7. **Do not** add READMEs, demos, or extra docs unless explicitly requested.

## Copyright Header

C++ headers and sources use a **Doxygen file block** at the top:

```cpp
/**
 * @file filename.hpp
 * @brief One-line description.
 * @copyright Copyright (c) YYYY Boyle Development Team
 */
```

Python (when used):

```python
# Copyright (c) YYYY Boyle Development Team. All rights reserved.
```

CMake scripts and `.mdc` rule files do not use this header where YAML frontmatter or CMake policy must be first line.

## TODO Format

```
// TODO(github_username): description [optional issue or PR reference]
```

Use `//` in C++ and `#` in Python. Include a real GitHub username in parentheses.

## File Creation Policy

- **Do not** create demo programs, sample apps, unsolicited README files, or documentation files unless the user explicitly asks.
- **Do not** add dependencies without using the project’s workflow (**CPM** / **vcpkg** for C++; **uv** for Python when enabled).
- Keep changes focused on the requested task.

## Security

- **Never** hardcode API keys, tokens, or passwords in source or committed config.
- **Never** commit `.env` files with real credentials.
- Prefer environment variables or secret managers for sensitive values.

---

# andrej-karpathy-skills

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:

- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:

- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:

- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:

- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.
