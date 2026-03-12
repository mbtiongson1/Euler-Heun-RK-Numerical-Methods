# Euler, Heun, Runge-Kutta Numerical Methods with Polynomial Fitting

Complete numerical analysis framework for solving first-order initial value problems (IVPs), systems of IVPs, and approximating solutions with polynomial fitting methods.

## Overview

This project implements:

**Numerical Integration Methods (Single IVP):**
- Euler Method
- Predictor-Corrector (Heun)
- Ralston (RK2.2)
- Runge-Kutta 3rd Order (RK3)
- Classical Runge-Kutta (RK4)

**Numerical Integration Methods (System IVP):**
- Heun Method for systems
- RK4 Method for systems
- Multi-method system comparison runner

**Polynomial Fitting Methods:**
- Vandermonde Matrix Method
- Lagrange Interpolation
- Least-Squares fitting comparison across selected methods

**Supporting Features:**
- Finite Difference derivative calculations
- Smart sampling and validation controls
- Optional actual/exact solution comparison
- CSV export of results
- Comparative visualization and analysis

---

## Project Structure

```
├── numerical_methods/              # Canonical Python package (import from here)
│   ├── problems/                   # Edit configs here (single IVP, systems, shooting, etc.)
│   ├── methods/                    # Individual method scripts (Euler/Heun/RK...)
│   ├── main/                       # Main drivers (solver/systems/shooting/function)
│   ├── fitting/                    # Vandermonde/Lagrange/least-squares modules
│   ├── fd/                         # Finite-difference workflows
│   ├── experiments/                # Scratch/analysis scripts
│   ├── utils.py
│   ├── matrix.py
│   └── paths.py                    # out/ path helpers
├── docs/
│   └── index.html
├── out/                            # generated outputs (ignored by git)
├── root/                           # optional wrapper scripts (python root/solver.py, etc.)
└── README.md
```

---

## Quick Start

### 1. Configure Your Problem

- Use `numerical_methods/problems/ivp.py` for single ODEs: `y' = f(x, y)`
- Use `numerical_methods/problems/ivpsystems.py` for systems:
  - `x' = f(x, y, t)`
  - `y' = g(x, y, t)`

### 2. Select Method(s)

- In `numerical_methods/problems/ivp.py`, set:
  - `method = 'euler' | 'heun' | 'rk22' | 'rk3' | 'rk4'`
  - `ls_methods = [...]` for least-squares comparisons
- In `numerical_methods/problems/ivpsystems.py`, set:
  - `method = 'heun' | 'rk4'` (for system-focused scripts)
- In `numerical_methods/main/systems.py`, use:
  - `active_methods = ['euler', 'heun', 'ral', 'rk3', 'rk4']` (any subset)

### 3. Run a Driver

```bash
# Preferred (module runs)
python -m numerical_methods.main.function
python -m numerical_methods.main.solver
python -m numerical_methods.main.systems
python -m numerical_methods.fd.FD

# Optional wrapper scripts (if you prefer python <file>.py style)
python root/function.py
python root/solver.py
python root/systems.py
python root/FD.py
```

### 4. Run Individual Methods

```bash
python -m numerical_methods.methods.euler
python -m numerical_methods.methods.heun
python -m numerical_methods.methods.rk22
python -m numerical_methods.methods.rk3
python -m numerical_methods.methods.rk4
python -m numerical_methods.methods.heunsystems
python -m numerical_methods.methods.rk4systems

# Optional wrappers
python root/rk4.py
```

Outputs that are written to disk (CSVs) are written under `out/csv/`.

---

## Numerical Methods

### Core Integration Methods (Single IVP)

All single-IVP methods solve `y' = f(x, y)` with initial condition `y(x0) = y0`.

**Implementation note:** counter-based loops are used to avoid floating-point drift in endpoint stepping.

#### Euler Method

```
y_{n+1} = y_n + h*f(x_n, y_n)
```

#### Heun Method (Predictor-Corrector)

```
Predictor:  y_p = y_n + h*f(x_n, y_n)
Corrector:  y_{n+1} = y_n + (h/2)*[f(x_n, y_n) + f(x_{n+1}, y_p)]
```

#### RK2.2 (Ralston)

```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + (3/4)h, y_n + (3/4)k1)
y_{n+1} = y_n + (1/3)k1 + (2/3)k2
```

#### RK3 (Kutta 3rd order)

```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h, y_n - k1 + 2k2)
y_{n+1} = y_n + (1/6)*(k1 + 4k2 + k3)
```

#### RK4 (Classical)

```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h/2, y_n + k2/2)
k4 = h*f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)*(k1 + 2k2 + 2k3 + k4)
```

### System IVP Methods

System scripts solve:

```
x' = f(x, y, t)
y' = g(x, y, t)
```

Implemented system solvers:
- Heun predictor-corrector (`heunsystems.py`)
- RK4 (`rk4systems.py`)
- Combined comparison with tables and plots (`numerical_methods/main/systems.py`)

---

## Polynomial Fitting Methods

### Data Sampling Strategy

To avoid overfitting and enforce valid polynomial selection:

1. Validate spacing/selection constraints
2. Sample points consistently across the interval
3. Build polynomial using selected points

### Vandermonde Matrix Method

Constructs and solves:

```
V * a = y
```

where `a` contains polynomial coefficients.

**Run standalone:**

```bash
python vandermonde.py
python vandermonde_manual.py
```

### Lagrange Interpolation

Constructs interpolation polynomial directly from selected data points.

**Run standalone:**

```bash
python lagrange.py
```

### Least-Squares Comparison

Fits degree-`p` polynomial models for methods listed in `ls_methods`.

**Run via:**

```bash
python function.py
```

---

## Finite Difference Calculations

Module `FD.py` computes finite-difference approximations from numerical method outputs.

Implemented formulas include:
- Forward difference
- Central difference
- Backward difference

---

## Configuration

### Single IVP (`numerical_methods/problems/ivp.py`)

```python
import math

def f(x, y):
    return math.sqrt(x) * math.sin(2*x) - 5*y

x0 = 0
y0 = 0
xn = 2.4

# choose one style
m = 1
h = 0.3 / (2 ** m)
# n = 9
# h = xn / (n - 1)

method = 'heun'                      # euler, heun, rk22, rk3, rk4
p = 10
ls_methods = ['euler', 'heun']

y_actual = None                      # or list of values, or generated values
```

### System IVP (`numerical_methods/problems/ivpsystems.py`)

```python
import math

def f(x, y, t):
    return x*y + t

def g(x, y, t):
    return y*t + x

x0 = 1
y0 = -1
t0 = 0
tn = 0.2

m = 0
h = 0.05 / (2 ** m)
# n = 9
# h = tn / (n - 1)

method = 'rk4'                       # heun or rk4
p = 10
ls_methods = ['heun', 'rk4']

x_actual, y_actual = None, None      # or generated/manual lists
```

---

## Main Orchestrators

### `numerical_methods/main/function.py`

Single-IVP end-to-end workflow:
1. Compute numerical solution for `method`
2. Build Vandermonde and Lagrange polynomials
3. Run least-squares across `ls_methods`
4. Plot and export results to `out/csv/output_fit.csv`

### `numerical_methods/main/solver.py`

Compares single-IVP methods in one run and can compute errors when `y_actual` is provided.

### `numerical_methods/main/systems.py`

Compares selected system-IVP methods (`active_methods`) with optional actual-solution overlays and error plots.

### `numerical_methods/main/shooting.py`

Compares selected shooting-problem methods (`active_methods`) with optional actual-solution overlays and error plots.

---

## Testing Framework (`test.py`)

Used to compare selected outputs across refinements and inspect consistency.

Run with:

```bash
python test.py
```

---

## CSV Export

Generated files may include:
- `out/csv/output.csv` (single-IVP multi-method comparison)
- `out/csv/output_fit.csv` (fitting workflow export)
- `out/csv/output_systems_x.csv` and `out/csv/output_systems_y.csv` (system-IVP comparison tables)
- `out/csv/output_systems_z.csv` (shooting workflow table, if produced)
- `out/csv/FD.csv` (finite difference table)

---

## Utilities (`utils.py`)

Helper functions include:
- `print_table(...)`
- `print_table_csv(...)`
- `plot_polynomial(...)`
- `plot_polynomials_compare(...)`

Utilities live at `numerical_methods/utils.py` and are imported as `numerical_methods.utils`.

---

## Example Usage

### Single IVP Example

1. Set in `numerical_methods/problems/ivp.py`:

```python
def f(x, y):
    return -2*x*y

x0 = 0
y0 = 1
xn = 2
method = 'rk4'
y_actual = None
```

2. Run:

```bash
python -m numerical_methods.main.solver
```

### System IVP Example

1. Set in `numerical_methods/problems/ivpsystems.py`:

```python
def f(x, y, t):
    return x + y

def g(x, y, t):
    return x - y

x0 = 1
y0 = 0
t0 = 0
tn = 1
```

2. Run:

```bash
python -m numerical_methods.main.systems
```

---

## Key Features

- Accurate step progression with counter-based loops
- Multiple numerical methods for single and system IVPs
- Optional exact/actual-solution error analysis
- Polynomial approximation (Vandermonde, Lagrange, least-squares)
- CSV + console + plot outputs for analysis
- Scripts usable as standalone runs for coursework workflows

---

## Troubleshooting

**`TypeError: object of type 'NoneType' has no len()`**
- Cause: Script expects actual data but `y_actual` (or `x_actual, y_actual`) is `None`.
- Fix: Use updated scripts and keep actual arrays as either valid lists or `None` consistently.

**Actual/reference arrays do not match expected length**
- Cause: Provided actual values do not align with computed checkpoints.
- Fix: Regenerate from helper functions using the same interval/step settings, or supply lists with consistent lengths.

**`Unknown method` error in `function.py`**
- Cause: `method` not in supported set for that script.
- Fix: Use one of `euler`, `heun`, `rk22`, `rk4` for `function.py`.

**Plots not appearing**
- Cause: Missing matplotlib or non-interactive backend.
- Fix: `pip install matplotlib` and run in a graphical session.

**Numerical instability or oscillation**
- Cause: Step size too large for the problem dynamics.
- Fix: Reduce `h` (or increase refinement `m`/`n`) and verify ODE/system definitions.

---

## Dependencies

- Python 3.7+
- NumPy
- Matplotlib
- `math` (standard library)

Install with:

```bash
pip install numpy matplotlib
```

---

## License

Academic use for numerical methods coursework.
