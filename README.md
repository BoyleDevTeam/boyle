<h1 align="center">Boyle</h1>

<p align="center">
  <b>The Fundamental Math Basis for Autonomous Driving Vehicles and Robotics.</b>
</p>

<p align="center">
  <a href="README.zh-cn.md">中文版</a> •
  <a href="LICENSE">BSD-3-Clause License</a> •
  <a href="https://github.com/BoyleDevTeam/boyle">GitHub</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-23-blue.svg" alt="C++ Standard">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg" alt="License">
</p>

---

## Introduction

Boyle is a high-performance C++23 mathematical library purpose-built for **autonomous driving** and **robotics** applications. In these domains, engineers routinely need:

- **Spline-based path representation** — smoothly interpolating waypoints into drivable trajectories with curvature, tangent, and arc-length queries.
- **Convex optimization** — solving quadratic programs in real-time for trajectory smoothing, obstacle avoidance, and control allocation.
- **Motion planning models** — combining reference lines, time–arc-length mappings, and hard/soft constraints (borders, fences) to generate feasible, comfortable motions.
- **Efficient linear algebra** — fast vector/matrix operations (with optional BLAS/LAPACK and AVX-512 SIMD) as the computational foundation for all of the above.

Boyle provides all of these capabilities in a unified, header-mostly library with modern C++23 concepts and templates. Each component is a fine-grained CMake target (e.g., `Boyle::math_vec2`, `Boyle::cvxopm_osqp_solver`) so you can link only what you need.

## Module Overview

```
src/boyle/
├── common/                         # boyle::common — Utilities
│   ├── fsm.hpp
│   └── utils/
│       ├── aligned_allocator.hpp
│       ├── aligned_memory_resource.hpp
│       ├── chrono_inspector.hpp
│       ├── exec_on_exit.hpp
│       ├── in_in_in_out_result.hpp
│       ├── logging.hpp
│       └── macros.hpp
├── math/                           # boyle::math — Mathematical Primitives
│   ├── concepts.hpp
│   ├── type_traits.hpp
│   ├── utils.hpp
│   ├── cubic_interpolation.hpp
│   ├── quintic_interpolation.hpp
│   ├── fft.hpp
│   ├── chebyshev.hpp
│   ├── dense/                      #   Dense Linear Algebra
│   │   ├── vec2.hpp
│   │   ├── vec3.hpp
│   │   ├── vector.hpp
│   │   ├── vectorx.hpp
│   │   ├── vector_view.hpp
│   │   ├── matrix.hpp
│   │   ├── matrixx.hpp
│   │   ├── matrix_view.hpp
│   │   ├── lu_dcmp.hpp
│   │   ├── qr_dcmp.hpp
│   │   ├── partial_pivot.hpp
│   │   ├── dense_traits.hpp
│   │   └── detail/
│   ├── sparse/                     #   Sparse Matrices
│   │   ├── dok_matrix.hpp
│   │   ├── lil_matrix.hpp
│   │   ├── coo_matrix.hpp
│   │   ├── csr_matrix.hpp
│   │   ├── csc_matrix.hpp
│   │   ├── sparse_matrix_proxy.hpp
│   │   ├── index_pair.hpp
│   │   └── sparse_traits.hpp
│   ├── curves/                     #   Piecewise Curves
│   │   ├── piecewise_linear_curve.hpp
│   │   ├── piecewise_cubic_curve.hpp
│   │   ├── piecewise_quintic_curve.hpp
│   │   ├── curve_proxy.hpp
│   │   ├── sl.hpp
│   │   └── slv.hpp
│   ├── functions/                  #   Piecewise Functions
│   │   ├── piecewise_linear_function.hpp
│   │   ├── piecewise_cubic_function.hpp
│   │   ├── piecewise_quintic_function.hpp
│   │   ├── function_proxy.hpp
│   │   └── boundary_mode.hpp
│   └── mdfunctions/                #   Multi-dimensional Functions
│       ├── linear_mdfunction.hpp
│       ├── quadratic_mdfunction.hpp
│       ├── rosenbrock_function.hpp
│       └── mdfunction_proxy.hpp
├── cvxopm/                         # boyle::cvxopm — Convex Optimization
│   ├── settings.hpp
│   ├── result.hpp
│   ├── info.hpp
│   ├── problems/
│   │   ├── qp_problem.hpp
│   │   └── dense_problem.hpp
│   └── solvers/
│       ├── osqp_solver.hpp / .cpp
│       ├── lbfgs_solver.hpp
│       ├── bfgs_solver.hpp
│       ├── amoeba_solver.hpp
│       ├── lnsrch_solver.hpp
│       └── detail/
└── kinetics/                       # boyle::kinetics — Motion Planning
    ├── route_line.hpp
    ├── path.hpp
    ├── motion.hpp
    ├── fence.hpp
    ├── border.hpp
    ├── dualism.hpp
    └── models/
        ├── route_line_cubic_acc_model.hpp
        ├── route_line_quintic_acc_model.hpp
        ├── route_line_cubic_offset_model.hpp
        └── route_line_quintic_offset_model.hpp
```

### `boyle::math` — Mathematical Primitives

#### Dense Linear Algebra (`math/dense/`)

| File | Target | Description |
|------|--------|-------------|
| `vec2.hpp` | `Boyle::math_vec2` | Fixed-size 2D vector `Vec2<T>` with `euclidean()`, `angle()`, `rotate()`, `crossProj()`, `normalized()`, optional AVX-512 SIMD |
| `vec3.hpp` | `Boyle::math_vec3` | Fixed-size 3D vector `Vec3<T>` with `cross()`, `crossProj()`, optional AVX-512 SIMD |
| `vector.hpp` | `Boyle::math_vector` | Fixed-size dense vector `Vector<T, N>` with `dot()`, `euclidean()`, `conjugated()`, optional BLAS/LAPACK |
| `vectorx.hpp` | `Boyle::math_vectorx` | Dynamic-size dense vector `VectorX<T>` (heap-allocated, PMR-aware) |
| `vector_view.hpp` | `Boyle::math_vector_view` | Non-owning view `VectorView<T>` over contiguous memory with stride support |
| `matrix.hpp` | `Boyle::math_matrix` | Fixed-size dense matrix `Matrix<T, NRows, NCols>` with row/column-major layout, optional BLAS/LAPACK |
| `matrixx.hpp` | `Boyle::math_matrixx` | Dynamic-size dense matrix `MatrixX<T>` (heap-allocated, PMR-aware) |
| `matrix_view.hpp` | `Boyle::math_matrix_view` | Non-owning view `MatrixView<T>` over matrix memory |
| `lu_dcmp.hpp` | `Boyle::math_lu_dcmp` | LU decomposition `LuDcmp<T>` with partial pivoting for solving linear systems |
| `qr_dcmp.hpp` | `Boyle::math_qr_dcmp` | QR decomposition `QrDcmp<T>` for least-squares and linear system solving |

#### Sparse Matrices (`math/sparse/`)

| File | Target | Description |
|------|--------|-------------|
| `dok_matrix.hpp` | `Boyle::math_dok_matrix` | Dictionary-of-keys format `DokMatrix<T, I>` for incremental construction |
| `lil_matrix.hpp` | `Boyle::math_lil_matrix` | List-of-lists format `LilMatrix<T, I>` for row-wise construction |
| `coo_matrix.hpp` | `Boyle::math_coo_matrix` | Coordinate format `CooMatrix<T, I>` (row/col/value triplets) |
| `csr_matrix.hpp` | `Boyle::math_csr_matrix` | Compressed sparse row `CsrMatrix<T, I>` for efficient row access |
| `csc_matrix.hpp` | `Boyle::math_csc_matrix` | Compressed sparse column `CscMatrix<T, I>` for efficient column access |
| `sparse_matrix_proxy.hpp` | `Boyle::math_sparse_matrix_proxy` | Type-erased proxy for sparse matrix polymorphism (via Microsoft Proxy) |

#### Piecewise Curves (`math/curves/`)

| File | Target | Description |
|------|--------|-------------|
| `piecewise_linear_curve.hpp` | `Boyle::math_piecewise_linear_curve` | Linear interpolation curve through `Vec2`/`Vec3` anchor points |
| `piecewise_cubic_curve.hpp` | `Boyle::math_piecewise_cubic_curve` | Cubic spline curve with `tangent()`, `normal()`, `curvature()`, `inverse()`, customizable boundary conditions |
| `piecewise_quintic_curve.hpp` | `Boyle::math_piecewise_quintic_curve` | Quintic spline curve with higher-order smoothness |
| `curve_proxy.hpp` | `Boyle::math_curve_proxy` | Type-erased proxy for curve types |
| `sl.hpp` / `slv.hpp` | (included in curve targets) | Arc-length–lateral-offset coordinate types `SlDuplet`, `SlvTriplet` |

#### Piecewise Functions (`math/functions/`)

| File | Target | Description |
|------|--------|-------------|
| `piecewise_linear_function.hpp` | `Boyle::math_piecewise_linear_function` | Piecewise linear scalar/vector function |
| `piecewise_cubic_function.hpp` | `Boyle::math_piecewise_cubic_function` | Cubic spline scalar function with `derivative()`, `integral()`, boundary conditions |
| `piecewise_quintic_function.hpp` | `Boyle::math_piecewise_quintic_function` | Quintic spline scalar function with higher-order smoothness |
| `function_proxy.hpp` | `Boyle::math_function_proxy` | Type-erased proxy for function types |
| `boundary_mode.hpp` | (included in cubic/quintic targets) | Boundary condition specification `BoundaryMode<T>` |

#### Multi-dimensional Functions (`math/mdfunctions/`)

| File | Target | Description |
|------|--------|-------------|
| `linear_mdfunction.hpp` | `Boyle::math_linear_mdfunction` | Linear multi-dimensional function `LinearMdFunction<T>` |
| `quadratic_mdfunction.hpp` | `Boyle::math_quadratic_mdfunction` | Quadratic multi-dimensional function `QuadraticMdFunction<T>` |
| `rosenbrock_function.hpp` | `Boyle::math_rosenbrock_function` | Rosenbrock benchmark function `RosenbrockFunction<T>` |
| `mdfunction_proxy.hpp` | `Boyle::math_mdfunction_proxy` | Type-erased proxy for multivariate functions with `eval()` and `gradient()` |

#### Other Math Utilities

| File | Target | Description |
|------|--------|-------------|
| `concepts.hpp` | `Boyle::math_concepts` | C++23 concepts: `ScalarArithmetic`, `VecArithmetic`, `MatArithmetic`, `Allocatory`, etc. |
| `type_traits.hpp` | `Boyle::math_type_traits` | Type traits for template metaprogramming |
| `utils.hpp` | `Boyle::math_utils` | Utilities: `lerp()`, `linspace()`, `hasDuplicates()`, `nearestUpperElement()` |
| `cubic_interpolation.hpp` | `Boyle::math_cubic_interpolation` | Cubic interpolation algorithm |
| `quintic_interpolation.hpp` | `Boyle::math_quintic_interpolation` | Quintic interpolation algorithm |
| `fft.hpp` | `Boyle::math_fft` | FFT via pocketfft (forward, backward, ortho normalization) |
| `chebyshev.hpp` | `Boyle::math_chebyshev` | Chebyshev polynomial approximation `Chebyshev<T>` |

---

### `boyle::cvxopm` — Convex Optimization

#### Problem Formulations (`cvxopm/problems/`)

| File | Target | Description |
|------|--------|-------------|
| `qp_problem.hpp` | `Boyle::cvxopm_qp_problem` | Sparse QP problem `QpProblem<T, I>`: minimize ½xᵀPx + qᵀx subject to l ≤ Ax ≤ u |
| `dense_problem.hpp` | `Boyle::cvxopm_dense_problem` | Unconstrained dense problem `DenseProblem<T>` wrapping an `MdFunctionProxy` objective |

#### Solvers (`cvxopm/solvers/`)

| File | Target | Description |
|------|--------|-------------|
| `osqp_solver.hpp` | `Boyle::cvxopm_osqp_solver` | QP solver `OsqpSolver<T, I>` via the OSQP library, supports warm start and polishing |
| `lbfgs_solver.hpp` | `Boyle::cvxopm_lbfgs_solver` | Limited-memory BFGS solver `LbfgsSolver<T>` for unconstrained optimization |
| `bfgs_solver.hpp` | `Boyle::cvxopm_bfgs_solver` | BFGS solver `BfgsSolver<T>` for unconstrained optimization |
| `amoeba_solver.hpp` | `Boyle::cvxopm_amoeba_solver` | Nelder–Mead (simplex) solver `AmoebaSolver<T>` for derivative-free optimization |
| `lnsrch_solver.hpp` | `Boyle::cvxopm_lnsrch_solver` | Line search solver `LnsrchSolver<T>` |

#### Supporting Types

| File | Target | Description |
|------|--------|-------------|
| `settings.hpp` | `Boyle::cvxopm_settings` | Solver settings `Settings<T, I>` |
| `result.hpp` | `Boyle::cvxopm_result` | Solution result `Result<T>` (primal/dual variables) |
| `info.hpp` | `Boyle::cvxopm_info` | Solver info `Info<T, I>` (iteration count, status) |

---

### `boyle::kinetics` — Motion Planning

| File | Target | Description |
|------|--------|-------------|
| `route_line.hpp` | `Boyle::kinetics_route_line` | Reference line `RouteLine<T>` as a piecewise quintic curve for autonomous driving |
| `path.hpp` | `Boyle::kinetics_path` | 2D path `Path<T>` with `tangent()`, `normal()`, `curvature()`, `inverse()` |
| `motion.hpp` | `Boyle::kinetics_motion` | Time–arc-length mapping `Motion<T>` with `velocity()`, `accel()`, `jerk()` |
| `fence.hpp` | `Boyle::kinetics_fence` | Time–space constraints `HardFence<T>` / `SoftFence<T>` for longitudinal planning |
| `border.hpp` | `Boyle::kinetics_border` | Spatial constraints `HardBorder<T>` / `SoftBorder<T>` with chirality for lateral planning |
| `dualism.hpp` | `Boyle::kinetics_dualism` | Enumerations: `Chirality` (LEFT/RIGHT), `Actio` (BLOCKING/YIELDING) |

#### Motion Planning Models (`kinetics/models/`)

| File | Target | Description |
|------|--------|-------------|
| `route_line_cubic_acc_model.hpp` | `Boyle::kinetics_route_line_cubic_acc_model` | Cubic acceleration model — generates longitudinal motion via QP optimization |
| `route_line_quintic_acc_model.hpp` | `Boyle::kinetics_route_line_quintic_acc_model` | Quintic acceleration model — smoother longitudinal motion via QP optimization |
| `route_line_cubic_offset_model.hpp` | `Boyle::kinetics_route_line_cubic_offset_model` | Cubic offset model — generates lateral path offset via QP optimization |
| `route_line_quintic_offset_model.hpp` | `Boyle::kinetics_route_line_quintic_offset_model` | Quintic offset model — smoother lateral path offset via QP optimization |

---

### `boyle::common` — Utilities

| File | Target | Description |
|------|--------|-------------|
| `fsm.hpp` | `Boyle::common_fsm` | Finite state machine with event dispatching, state stack, and `FSM_INITIAL_STATE` macro |
| `logging.hpp` | `Boyle::common_logging` | spdlog-based logging macros (`BOYLE_LOG_TRACE`, `BOYLE_LOG_DEBUG`, …) |
| `macros.hpp` | `Boyle::common_macros` | Convenience macros for copy/move/singleton patterns |
| `aligned_allocator.hpp` | `Boyle::common_aligned_allocator` | 32-byte aligned allocator `AlignedAllocator<T>` |
| `aligned_memory_resource.hpp` | `Boyle::common_aligned_memory_resource` | PMR aligned memory resource |
| `exec_on_exit.hpp` | `Boyle::common_exec_on_exit` | RAII cleanup handler `ExecOnExit` |
| `chrono_inspector.hpp` | `Boyle::common_chrono_inspector` | Execution time measurement |

---

## Prerequisites

The following tools must be installed on your system:

| Tool | Purpose |
|------|---------|
| **CMake** (≥ 3.31) | Build system generator |
| **Ninja** | Build backend (used by all presets) |
| **GCC** (≥ 14) or **Clang** (≥ 20) | C++23-capable compiler |
| **clangd** | Language server for IDE integration |
| **clang-tidy** | Static analysis |
| **clang-format** | Code formatting |
| **ccache** | Compilation cache for faster rebuilds |

All other dependencies (Boost, OSQP, spdlog, doctest, etc.) are automatically managed by [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) and do not require manual installation.

## Build

```bash
# Configure (Linux + GCC)
cmake --preset linux-gcc-x64

# Build
cmake --build --preset linux-gcc-x64

# Run tests
ctest --preset linux-gcc-x64
```

For Clang:

```bash
cmake --preset linux-clang-x64
cmake --build --preset linux-clang-x64
ctest --preset linux-clang-x64
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `BOYLE_CHECK_PARAMS` | `OFF` | Enable runtime parameter validation |
| `BOYLE_BUILD_TESTING` | `ON` | Build unit tests |
| `BOYLE_ENABLE_INSTALL` | `ON` | Enable install targets |
| `BOYLE_USE_SIMD` | `OFF` | SIMD implementation (`AVX512` or `OFF`) |
| `BOYLE_USE_BLAS_LAPACK` | `OpenBLAS` | BLAS/LAPACK backend (`OpenBLAS`, `Netlib`, `MKL`, `OFF`) |
| `BOYLE_USE_BOOST_UNORDERED` | `ON` | Use Boost.Unordered containers |
| `BUILD_SHARED_LIBS` | `OFF` | Build shared libraries |

### Install

```bash
cmake --install out/build/linux-gcc-x64 --prefix /your/install/path
```

After installation, use Boyle in your own CMake project by linking the fine-grained targets you need:

```cmake
find_package(Boyle REQUIRED)
target_link_libraries(your_target
  PRIVATE
    Boyle::math_vec2
    Boyle::math_piecewise_cubic_curve
    Boyle::cvxopm_osqp_solver
    Boyle::kinetics_route_line
)
```

## Usage Examples

### 2D Vector Operations

```cpp
#include "boyle/math/dense/vec2.hpp"

using namespace boyle::math;

// Construct a 2D vector
Vec2d a(1.0, std::numbers::sqrt3);

// Basic queries
double len   = a.euclidean();             // 2.0
double angle = a.angle();                 // π/3
Vec2d  unit  = a.normalized();            // (0.5, √3/2)

// Rotation
Vec2d rot = a.rotate(M_PI / 6.0);        // (0, 2)

// Cross product projection & dot product
double cross = a.crossProj(Vec2d{0.0, 1.0});
double dot   = a.dot(Vec2d{1.0, 0.0});   // 1.0
```

### Fixed-size and Dynamic Vectors

```cpp
#include "boyle/math/dense/vector.hpp"
#include "boyle/math/dense/vectorx.hpp"

using namespace boyle::math;

// Fixed-size vector
Vector<double, 4> v;
v[0] = 1.0; v[1] = 2.0; v[2] = 3.0; v[3] = 4.0;

double norm = v.euclidean();
double d    = v.dot(v);
bool   same = v.identicalTo(v * 1.000000001);  // true (within tolerance)
auto   conj = v.conjugated();

// Dynamic-size vector (PMR-aware)
VectorX<double> dyn(100);
```

### Cubic Spline Curves

```cpp
#include "boyle/math/curves/piecewise_cubic_curve.hpp"
#include "boyle/math/dense/vec2.hpp"
#include "boyle/math/functions/boundary_mode.hpp"

using namespace boyle::math;

// Create anchor points on a semicircle (radius = 2)
std::vector<Vec2d> anchors(21);
for (std::size_t i = 0; i < 21; ++i) {
    double theta = M_PI * i / 20.0;
    anchors[i] = Vec2d{2.0 * std::cos(theta), 2.0 * std::sin(theta)};
}

// Set clamped boundary conditions (first-order derivative)
BoundaryMode<Vec2d> b0{.order = 1, .derivative = Vec2d{0.0, 1.0}};
BoundaryMode<Vec2d> bf{.order = 1, .derivative = Vec2d{0.0, -1.0}};

// Construct the cubic spline curve
PiecewiseCubicCurve<Vec2d> curve{anchors, b0, bf};

// Query the curve at arc-length s = 1.0
Vec2d  point   = curve(1.0);
Vec2d  tangent = curve.tangent(1.0);
Vec2d  normal  = curve.normal(1.0);
double kappa   = curve.curvature(1.0);    // ≈ 0.5 (1/radius)

// Cartesian → Frenet (s, l) coordinate conversion
SlDupletd sl        = curve.inverse(Vec2d{1.0, 1.0});
Vec2d     recovered = curve(sl);           // (s, l) → Cartesian
```

### Cubic Spline Functions

```cpp
#include "boyle/math/functions/piecewise_cubic_function.hpp"
#include "boyle/math/functions/boundary_mode.hpp"
#include "boyle/math/utils.hpp"

using namespace boyle::math;

// Generate sample points for a polynomial
std::vector<double> ts = linspace(-2.0, 2.0, 41);
std::vector<double> ys;
for (double t : ts) {
    ys.push_back(0.45 + 5.3*t - 1.3*t*t + 0.65*t*t*t);
}

// Boundary conditions: natural (second-order derivative)
BoundaryMode<double> b0{.order = 2, .derivative = -2.6};
BoundaryMode<double> bf{.order = 2, .derivative =  1.3};

PiecewiseCubicFunction<double> func{ts, ys, b0, bf};

double val  = func(0.5);                  // f(0.5)
double dval = func.derivative(0.5);       // f'(0.5)
double area = func.integral(-1.0, 1.0);   // ∫₋₁¹ f(t) dt
```

### Quadratic Programming (OSQP)

```cpp
#include "boyle/cvxopm/problems/qp_problem.hpp"
#include "boyle/cvxopm/solvers/osqp_solver.hpp"
#include "boyle/cvxopm/settings.hpp"

using namespace boyle::cvxopm;

// Minimize: 2x₀² + x₁² + x₀x₁ + x₀ + x₁
// Subject to: x₀ + x₁ = 1, 0 ≤ x₀ ≤ 0.7, 0 ≤ x₁ ≤ 0.7
QpProblem<double, int> qp(2, 3);   // 2 variables, 3 constraints
qp.updateQuadCostTerm(0, 0, 2.0);
qp.updateQuadCostTerm(1, 1, 1.0);
qp.updateQuadCostTerm(0, 1, 1.0);
qp.updateLinCostTerm(0, 1.0);
qp.updateLinCostTerm(1, 1.0);
qp.updateConstrainTerm(0, {{0, 1.0}, {1, 1.0}}, 1.0, 1.0);
qp.updateConstrainTerm(1, {{0, 1.0}}, 0.0, 0.7);
qp.updateConstrainTerm(2, {{1, 1.0}}, 0.0, 0.7);

OsqpSolver<double, int> solver{Settings<double, int>{.polishing = true}};
auto [result, info] = solver.solve(qp);
// result.prim_vars[0] ≈ 0.3, result.prim_vars[1] ≈ 0.7

// Warm start with a known initial point
std::vector<double> prim_vars_0{0.3, 0.7};
std::vector<double> dual_vars_0{-2.9, 0.0, 0.2};
auto [result2, info2] = solver.solve(qp, prim_vars_0, dual_vars_0);
```

### Unconstrained Optimization (L-BFGS)

```cpp
#include "boyle/cvxopm/solvers/lbfgs_solver.hpp"
#include "boyle/cvxopm/problems/dense_problem.hpp"
#include "boyle/math/mdfunctions/rosenbrock_function.hpp"
#include "boyle/math/dense/vectorx.hpp"

using namespace boyle::cvxopm;
using namespace boyle::math;

// Minimize the Rosenbrock function: f(x,y) = (a-x)² + b(y-x²)²
DenseProblem<double> problem(
    makeMdFunctionProxy(RosenbrockFunction<pmr::VectorX<double>>(1.0, 100.0))
);

LbfgsSolver<double> solver{.settings{.eps_abs{1E-6}}};
pmr::VectorX<double> x0({-1.0, 1.0});
auto [result, info] = solver.solve(problem, std::move(x0));
// result.prim_vars ≈ (1.0, 1.0), converges within ~40 iterations
```

### Sparse Matrices

```cpp
#include "boyle/math/sparse/dok_matrix.hpp"
#include "boyle/math/sparse/csr_matrix.hpp"

using namespace boyle::math;

// Build a sparse matrix incrementally with DOK format
DokMatrix<double, int> dok(3, 3);
dok.updateCoeff(0, 0, 1.0);
dok.updateCoeff(1, 1, 2.0);
dok.updateCoeff(2, 2, 3.0);
dok.updateCoeff(0, 2, 0.5);

// Convert to CSR for efficient row-wise access
CsrMatrix<double, int> csr(dok);
```

## Namespaces

| Namespace | Description |
|-----------|-------------|
| `boyle::math` | Mathematical primitives: vectors, matrices, curves, functions, FFT, Chebyshev |
| `boyle::math::detail` | Internal implementation details |
| `boyle::cvxopm` | Convex optimization: QP problems, OSQP / L-BFGS / BFGS / Amoeba solvers |
| `boyle::kinetics` | Motion planning: paths, route lines, motions, borders, fences, planning models |
| `boyle::common` | Utilities: FSM, logging, aligned memory, macros |

## License

This project is licensed under the [BSD 3-Clause License](LICENSE).

## Acknowledgments

Boyle is built upon the following outstanding open-source projects. We extend our sincere gratitude to their developers and communities:

| Project | Version | Description | Homepage |
|---------|---------|-------------|----------|
| [Boost](https://www.boost.org/) | 1.90.0 | Serialization, unordered containers, hash | https://www.boost.org/ |
| [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) | 0.42.1 | CMake package manager | https://github.com/cpm-cmake/CPM.cmake |
| [cxxopts](https://github.com/jarro2783/cxxopts) | 3.3.1 | Command-line option parsing | https://github.com/jarro2783/cxxopts |
| [doctest](https://github.com/doctest/doctest) | 2.4.12 | Fast C++ testing framework | https://github.com/doctest/doctest |
| [Matplot++](https://github.com/alandefreitas/matplotplusplus) | 1.2.2 | C++ graphics library for data visualization | https://github.com/alandefreitas/matplotplusplus |
| [Microsoft Proxy](https://github.com/microsoft/proxy) | 4.0.2 | Type-erased polymorphism library | https://github.com/microsoft/proxy |
| [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) | 0.3.32 | Optimized BLAS/LAPACK implementation (optional) | https://github.com/OpenMathLib/OpenBLAS |
| [OSQP](https://osqp.org/) | 1.0.0 | Quadratic programming solver | https://osqp.org/ |
| [pocketfft](https://github.com/mreineck/pocketfft) | cpp branch | Header-only FFT library | https://github.com/mreineck/pocketfft |
| [qdldl](https://github.com/osqp/qdldl) | 0.1.9 | Sparse LDL factorization (OSQP dependency) | https://github.com/osqp/qdldl |
| [spdlog](https://github.com/gabime/spdlog) | 1.17.0 | Fast C++ logging library | https://github.com/gabime/spdlog |
