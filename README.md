# Euler, Heun, Runge-Kutta Numerical Methods with Polynomial Fitting

Complete numerical analysis framework for solving first-order initial value problems (IVPs) and approximating solutions with polynomial fitting methods.

## Overview

This project implements:

**Numerical Integration Methods:**
- Euler Method
- Predictor–Corrector (Heun)
- Ralston (RK2.2)
- Classical Runge-Kutta (RK4)

**Polynomial Fitting Methods:**
- Vandermonde Matrix Method (least-squares fitting)
- Lagrange Interpolation (exact interpolation)

**Supporting Features:**
- Finite Difference Derivative Calculations (dy/dx, d2y/dx2)
- Smart Data Point Sampling & Validation
- CSV Export of Results
- Comparative Visualization & Analysis

---

## Project Structure

```
├── ivp.py                  # Central configuration file
├── euler.py                # Euler method implementation
├── heun.py                 # Heun method implementation
├── rk22.py                 # RK2 (Ralston) implementation
├── rk4.py                  # RK4 method implementation
├── FD.py                   # Finite difference calculations
├── vandermonde.py          # Vandermonde polynomial fitting
├── lagrange.py             # Lagrange interpolation
├── function.py             # Main orchestrator with validation
├── utils.py                # Utility functions (plotting, display)
├── test.py                 # Comparison testing framework
├── output.csv              # CSV output (generated)
└── README.md               # This file
```

---

## Quick Start

### 1. Configure Your Problem (ivp.py)

Edit `ivp.py` to set your initial value problem:

```python
def f(x, y):
    return math.cos(x) - x*math.sin(y)  # Your differential equation

x0 = 0                  # Initial x
y0 = 1                  # Initial y
xn = 1.5                # End of interval
n = 9                   # Number of data points for n-refinement
p = 4                   # Polynomial degree for fitting
method = 'rk4'          # Numerical method: euler, heun, rk22, rk4
```

**Important:** Ensure `(n-1)` is a multiple of `p`:
- Example: n=9, p=4 → (9-1)=8, and 8 % 4 = 0 ✓
- If invalid, the system will suggest valid combinations

### 2. Run the Main Solver

```bash
# Run complete analysis with polynomial fitting
python function.py

# Run individual methods
python rk4.py              # Just RK4
python vandermonde.py      # Vandermonde fitting only
python lagrange.py         # Lagrange fitting only
python test.py             # Compare methods across refinement levels
```

### 3. View Results

- **Console Output:** Matrices, coefficients, polynomial equations
- **Plots:** Three matplotlib windows comparing methods
- **CSV Export:** `output.csv` with detailed data

---

## Numerical Methods

### Core Integration Methods

All methods solve `y' = f(x, y)` with initial condition `y(x0) = y0`.

**Implementation Note:** All methods use counter-based loops to prevent floating-point accumulation:
```python
for n in range(num_steps + 1):
    x = x0 + n * h  # Recalculate from origin, not accumulation
```

This ensures exact termination at `xn` without overshoot.

#### Euler Method
Simplest method, O(h) error per step.
```
y_{n+1} = y_n + h*f(x_n, y_n)
```

#### Heun Method (Predictor-Corrector)
Improved Euler with two evaluations, O(h²) error per step.
```
Predictor:  y_p = y_n + h*f(x_n, y_n)
Corrector:  y_{n+1} = y_n + (h/2)*[f(x_n, y_n) + f(x_{n+1}, y_p)]
```

#### RK2 (Ralston Method)
Alternative second-order method.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + (2/3)*h, y_n + (2/3)*k1)
y_{n+1} = y_n + (1/4)*k1 + (3/4)*k2
```

#### RK4 (Classical Runge-Kutta)
Most popular, O(h⁴) error per step.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h/2, y_n + k2/2)
k4 = h*f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
```

---

## Polynomial Fitting Methods

### Data Sampling Strategy

To avoid over-fitting and ensure proper polynomial degree usage:

1. **Validation:** Check if `(n-1) % p == 0`
2. **Sampling:** Select every k-th point where `k = (n-1) / p`
3. **Result:** Exactly `p+1` points spanning the full range, including endpoints

Example: n=9, p=4
- Intervals: 8 (from index 0 to 8)
- Step: 8/4 = 2
- Sampled indices: [0, 2, 4, 6, 8]
- Points: [x₀, x₂, x₄, x₆, x₈] covering the full range

### Vandermonde Matrix Method

Builds and solves a Vandermonde matrix for least-squares polynomial fitting.

**Matrix Form:**
```
V = [1   x₀   x₀²   ... x₀ᵖ]      a = [a₁]      y = [y₀]
    [1   x₁   x₁²   ... x₁ᵖ]          [a₂]          [y₁]
    [... ... ...  ... ... ]           [...]          [...] 
    [1   xₚ   xₚ²   ... xₚᵖ]          [aₚ₊₁]         [yₚ]

Solve: V·a = y using least-squares (NumPy lstsq)
```

**Output Format:**
```
Polynomial: f(x) = a₁ + a₂*x + a₃*x² + ... + aₚ₊₁*xᵖ
```

**Run Standalone:**
```bash
python vandermonde.py
```

### Lagrange Interpolation

Constructs a polynomial that passes exactly through all selected data points.

**Basis Functions:**
```
L_i(x) = ∏_{j≠i} (x - xⱼ)/(xᵢ - xⱼ)

P(x) = ∑ yᵢ * L_i(x)
```

**Advantages:**
- Exact interpolation (passes through all points)
- No matrix solving required
- Numerically stable for reasonable polynomial degrees

**Run Standalone:**
```bash
python lagrange.py
```

---

## Finite Difference Calculations

Module `FD.py` computes numerical derivatives using finite differences.

**Implemented Methods:**
- Forward difference: `f'(x) ≈ [f(x+h) - f(x)] / h`
- Central difference: `f'(x) ≈ [f(x+h) - f(x-h)] / (2h)`
- Backward difference: `f'(x) ≈ [f(x) - f(x-h)] / h`

**Usage:**
```python
from FD import compute_rk4
xs, ys = compute_rk4()
# xs: x values, ys: y values at those points
```

Each solver module (`euler.py`, `heun.py`, `rk22.py`, `rk4.py`) exports a `compute_*` function that returns `(xs, ys)` tuples.

---

## Configuration (ivp.py)

Central configuration file controlling all solvers and fitting methods.

### Parameters

```python
# Differential Equation
def f(x, y):
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0          # Starting x
y0 = 1          # Starting y
xn = 1.5        # Ending x

# Step Size (choose one method)
m = 0           # m-refinement: h = 0.3 / (2^m)
n = 9           # n-refinement: h = xn / (n-1)

# Numerical Method & Polynomial Degree
method = 'rk4'  # euler, heun, rk22, rk4
p = 4           # Polynomial degree (p+1 = number of fitting points)
```

### Validation Function

```python
validate_data_points(num_points, polynomial_degree)
```

- Checks if `(num_points - 1) % polynomial_degree == 0`
- Returns sampling indices for valid configurations
- Raises error with suggestions for invalid configurations

---

## Main Orchestrator (function.py)

Combines all methods and produces comprehensive output.

**Workflow:**
1. Compute numerical solution using selected method
2. Validate data points and sample appropriately
3. Build Vandermonde matrix
4. Solve using Vandermonde and Lagrange methods
5. Generate comparison plots
6. Export to CSV

**Output:**
```
Method: RK4, data points: 5, polynomial degree p=4

Vandermonde Matrix (all rows, coefficients only):
  [1.0000, 0.0000, 0.0000, 0.0000, 0.0000]  (x=0.0000, y=1.0000)
  [1.0000, 0.3750, 0.1406, 0.0527, 0.0198]  (x=0.3750, y=1.3005)
  ...

Coefficients:
a1 = 1.000000
a2 = 0.990219
a3 = -0.384869
a4 = -0.354567
a5 = 0.102217

Vandermonde Polynomial:
f(x) = 1.000000 + 0.990219*x^1 + -0.384869*x^2 + -0.354567*x^3 + 0.102217*x^4
```

**Run:**
```bash
python function.py
```

---

## Testing Framework (test.py)

Compare RK4 results across different refinement levels to verify consistency.

**Example: n=9 vs n=5**

```
RK4 COMPARISON: n=9 vs n=5
   n |        x |         y (n=9) |    n |         y (n=5)
------------------------------------------------------------
   0 |   0.0000 |        1.000000 |    0 |        1.000000
   1 |   0.1875 |        1.170618 |    - |               -
   2 |   0.3750 |        1.300533 |    1 |        1.300520
   ...
   8 |   1.5000 |        0.940182 |    4 |        0.940412

DIFFERENCES (where both have values)
       x |          y(n=9) |          y(n=5) |         |Δy| | Rel. Error %
----------------------------------------------------------------------
  0.0000 |        1.000000 |        1.000000 |     0.000000 |     0.000000
  0.3750 |        1.300533 |        1.300520 |     0.000013 |     0.001001
  0.7500 |        1.408935 |        1.408962 |     0.000028 |     0.001982
  1.1250 |        1.285786 |        1.285930 |     0.000145 |     0.011251
  1.5000 |        0.940182 |        0.940412 |     0.000230 |     0.024419
```

Shows alignment of common x values and error analysis.

**Run:**
```bash
python test.py
```

---

## CSV Export (output.csv)

Comprehensive data export with sections:

1. **Data Points:** All (x, y) values from numerical solver
2. **Vandermonde Matrix:** Full matrix with headers
3. **Vandermonde Coefficients:** a₁, a₂, ..., aₚ₊₁
4. **Lagrange Coefficients:** b₁, b₂, ..., bₚ₊₁

Import into spreadsheet software for further analysis.

---

## Utilities (utils.py)

Helper functions used throughout the project:

- `print_table(headers, rows)` - Pretty-print tabular data
- `plot_polynomial(xs, ys, coeffs, ...)` - Plot single polynomial fit
- `plot_polynomials_compare(xs, ys, coeff_lists, ...)` - Compare multiple fits

---

## Example Usage

### Solving a Custom ODE

1. Edit `ivp.py`:
```python
def f(x, y):
    return -2*x*y  # Your ODE: y' = -2xy

x0 = 0
y0 = 1
xn = 2.0
n = 9
p = 4
method = 'rk4'
```

2. Check validity: (9-1)=8, 8 % 4 = 0 ✓

3. Run:
```bash
python function.py
```

4. View results in console and plots

### Comparing Methods

```bash
# Show individual method outputs
python euler.py
python heun.py
python rk22.py
python rk4.py

# Compare fitting approaches
python vandermonde.py
python lagrange.py

# Test consistency
python test.py
```

---

## Key Features

✅ **Accurate Computation** - Counter-based loops eliminate floating-point drift  
✅ **Smart Sampling** - Validates (n-1) % p == 0 and samples appropriately  
✅ **Dual-Execution** - All modules work as importable libraries or standalone scripts  
✅ **Comprehensive Output** - Console display, plots, and CSV export  
✅ **Error Analysis** - Compares methods and shows relative differences  
✅ **Educational** - Clear documentation and test cases  

---

## Troubleshooting

**Error: "n % (p+1) is not zero"**
- Solution: Adjust n or p so (n-1) is divisible by p
- Example: n=9, p=4 works (8 % 4 = 0)

**Plots not appearing**
- Ensure matplotlib is installed: `pip install matplotlib`
- Use `plt.show()` or interactive notebook mode

**Numerical instability**
- Try reducing step size (increase n or m-refinement)
- Check that your ODE function is well-behaved

---

## Dependencies

- Python 3.7+
- NumPy (numerical operations)
- Matplotlib (plotting)
- Math (standard library)

Install with:
```bash
pip install numpy matplotlib
```

---

## License

Academic use for numerical methods coursework.

---

This framework provides a complete environment for exploring numerical methods and polynomial approximation in a unified, user-friendly package.
