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
├── ivp.py                  # Single-IVP configuration
├── ivpsystems.py           # System-IVP configuration
├── euler.py                # Euler method (single IVP)
├── heun.py                 # Heun method (single IVP)
├── rk22.py                 # RK2.2 / Ralston (single IVP)
├── rk3.py                  # RK3 (single IVP)
├── rk4.py                  # RK4 (single IVP)
├── heunsystems.py          # Heun method (system IVP)
├── rk4systems.py           # RK4 method (system IVP)
├── mainsolver.py           # Multi-method single-IVP comparison
├── mainsystems.py          # Multi-method system-IVP comparison
├── FD.py                   # Finite difference calculations
├── vandermonde.py          # Vandermonde polynomial fitting
├── vandermonde_manual.py   # Manual Vandermonde workflow
├── lagrange.py             # Lagrange interpolation
├── leastsquares.py         # Least-squares solver
├── function.py             # Combined fitting workflow for single IVP
├── utils.py                # Utility functions (tables, plotting)
├── test.py                 # Consistency/comparison testing
└── README.md               # This file
```

---

## Quick Start

### 1. Configure Your Problem

- Use `ivp.py` for single ODEs: `y' = f(x, y)`
- Use `ivpsystems.py` for systems:
  - `x' = f(x, y, t)`
  - `y' = g(x, y, t)`

### 2. Select Method(s)

- In `ivp.py`, set:
  - `method = 'euler' | 'heun' | 'rk22' | 'rk3' | 'rk4'`
  - `ls_methods = [...]` for least-squares comparisons
- In `ivpsystems.py`, set:
  - `method = 'heun' | 'rk4'` (for system-focused scripts)
- In `mainsystems.py`, use:
  - `active_methods = ['euler', 'heun', 'ral', 'rk3', 'rk4']` (any subset)

### 3. Run a Driver

```bash
# Single-IVP fitting workflow (Vandermonde, Lagrange, least-squares)
python function.py

# Single-IVP multi-method comparison table/plots
python mainsolver.py

# System-IVP multi-method comparison table/plots
python mainsystems.py

# Finite-difference workflow
python FD.py
```

### 4. Run Individual Methods

```bash
python euler.py
python heun.py
python rk22.py
python rk3.py
python rk4.py
python heunsystems.py
python rk4systems.py
```

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
- Combined comparison with tables and plots (`mainsystems.py`)

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

### Single IVP (`ivp.py`)

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

### System IVP (`ivpsystems.py`)

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

### `function.py`

Single-IVP end-to-end workflow:
1. Compute numerical solution for `method`
2. Build Vandermonde and Lagrange polynomials
3. Run least-squares across `ls_methods`
4. Plot and export results to `output.csv`

### `mainsolver.py`

Compares single-IVP methods in one run and can compute errors when `y_actual` is provided.

### `mainsystems.py`

Compares selected system-IVP methods (`active_methods`) with optional actual-solution overlays and error plots.

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
- `output.csv` (single-IVP fitting/comparison data)
- `output_systems_x.csv` and `output_systems_y.csv` (system-IVP comparison tables)

---

## Utilities (`utils.py`)

Helper functions include:
- `print_table(...)`
- `print_table_csv(...)`
- `plot_polynomial(...)`
- `plot_polynomials_compare(...)`

---

## Example Usage

### Single IVP Example

1. Set in `ivp.py`:

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
python mainsolver.py
```

### System IVP Example

1. Set in `ivpsystems.py`:

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
python mainsystems.py
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
