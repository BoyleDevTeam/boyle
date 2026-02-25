# Project Overview

Boyle is a **high-performance C++23** mathematical library for autonomous driving and robotics, distributed under the **BSD 3-Clause** license. It emphasizes dense and sparse linear algebra, piecewise curves and scalar/vector functions, convex optimization solvers, and kinetics-oriented motion models.

**Python bindings and tooling are planned** (packaging via **uv**, lint/format via **ruff**, tests via **pytest**). Until Python sources land in-tree, treat Python guidance as the target workflow for future modules.

## Module architecture

Dependency flow (strictly acyclic):

```
common → math → cvxopm → kinetics
```

- **common** — Shared utilities: finite-state helpers, allocators, logging (e.g. spdlog), macros (`BOYLE_*`), and low-level helpers.
- **math** — `dense`, `sparse`, `curves`, `functions`, `mdfunctions`, and related CMake/xmake targets.
- **cvxopm** — Convex optimization: problem formulations and solvers (OSQP, L-BFGS, etc., per dependencies).
- **kinetics** — Models used in motion planning contexts.

Implementation is predominantly **header-first** with fine-grained targets per class or small component.

## Dependencies

C++ dependencies are managed with **CPM** and **vcpkg** (vendored). Typical libraries include Boost (serialization, unordered), doctest, spdlog, OSQP/qdldl, OpenBLAS, and test helpers (cxxopts, Matplot++). Exact pins live under `cmake/` and `xmake/third_party/`.

## Library design principles

- **Namespaces** partition concepts: `boyle::math`, `boyle::cvxopm`, `boyle::kinetics`, `boyle::common`, with `detail` and `pmr` sub-namespaces where appropriate.
- **No circular** link dependencies between modules.
- **Exceptions** communicate misuse or numerical failure; optional **`BOYLE_CHECK_PARAMS`** adds inexpensive checks in debug-oriented builds.
- **Serialization**: Boost.Serialization where types opt in—always cover with roundtrip tests.

For command-line workflows and directory layout, see `02-build-and-test.md` and repository `AGENTS.md`.
