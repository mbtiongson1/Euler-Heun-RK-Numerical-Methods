# Numerical Differentiation via Finite Differences

Complete implementation of numerical differentiation methods using finite differences to approximate derivatives of solutions from numerical integration solvers.

## Overview

This branch focuses on **numerical differentiation** - computing approximate derivatives from discrete numerical solutions.

**Implemented Methods:**
- Euler Method
- Predictor–Corrector (Heun)
- Ralston (RK2.2)
- Classical Runge-Kutta (RK4)

**Derivative Approximations:**
- Forward Difference (at start)
- Central Difference (at interior points)
- Backward Difference (at end)

---

## Project Structure

```
├── ivp.py              # Central configuration file
├── euler.py            # Euler method implementation
├── heun.py             # Heun method implementation
├── rk22.py             # RK2 (Ralston) implementation
├── rk4.py              # RK4 method implementation
├── FD.py               # Finite difference derivative calculations
├── mainsolver.py       # Main orchestrator
├── utils.py            # Utility functions (plotting, display)
├── test.py             # Testing and analysis
├── output.csv          # CSV output (generated)
├── Figure 1 - Numerical Methods Comparison (y-values).png
├── Figure 2 - Percent Relative Errors Comparison.png
└── README.md           # This file
```

---

## Quick Start

### 1. Configure Your Problem (ivp.py)

Edit `ivp.py` to set your initial value problem:

```python
def f(x, y):
    """Your differential equation: y' = f(x, y)"""
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0                  # Starting x
y0 = 1                  # Starting y
xn = 1.5                # Ending x

# Step Size
h = 0.3                 # Fixed step size (or use refinement parameters)

# Exact Solution (optional)
def y_actual_func(x):
    return x*math.cos(x) + 1
```

### 2. Run the Main Solver

```bash
# Complete analysis with all methods
python mainsolver.py

# Individual methods
python euler.py
python heun.py
python rk22.py
python rk4.py

# Finite difference analysis
python FD.py

# Testing
python test.py
```

### 3. View Results

- **Console Output:** Tables with y-values, derivatives, and errors
- **Plots:** Comparison figures for y-values and percent relative errors
- **CSV Export:** `output.csv` with detailed numerical data

---

## Numerical Integration Methods

All methods solve `y' = f(x, y)` with initial condition `y(x0) = y0`.

### Euler Method
Simplest first-order method.
```
y_{n+1} = y_n + h*f(x_n, y_n)
```

### Heun Method (Predictor-Corrector)
Improved Euler with two function evaluations.
```
Predictor:  y_p = y_n + h*f(x_n, y_n)
Corrector:  y_{n+1} = y_n + (h/2)*[f(x_n, y_n) + f(x_{n+1}, y_p)]
```

### RK2 (Ralston Method)
Alternative second-order Runge-Kutta method.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + (2/3)*h, y_n + (2/3)*k1)
y_{n+1} = y_n + (1/4)*k1 + (3/4)*k2
```

### RK4 (Classical Runge-Kutta)
Most accurate, O(h⁴) local truncation error.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h/2, y_n + k2/2)
k4 = h*f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
```

---

## Numerical Differentiation

Once a solution `{(x_i, y_i)}` is obtained from a numerical solver, finite differences are used to approximate derivatives.

### Finite Difference Formulas

**Forward Difference (first point):**
```
f'(x_0) ≈ (y_1 - y_0) / h
O(h) error
```

**Central Difference (interior points):**
```
f'(x_i) ≈ (y_{i+1} - y_{i-1}) / (2h)
O(h²) error - most accurate
```

**Backward Difference (last point):**
```
f'(x_n) ≈ (y_n - y_{n-1}) / h
O(h) error
```

**Second Derivative (interior points):**
```
f''(x_i) ≈ (y_{i+1} - 2*y_i + y_{i-1}) / h²
O(h²) error
```

### Implementation

File `FD.py` provides compute functions for each method:

```python
from FD import compute_euler, compute_heun, compute_rk22, compute_rk4

xs, ys = compute_rk4()  # Returns (x_values, y_values)
```

Each solver automatically computes:
- **First Derivative (dy/dx)** at each point
- **Second Derivative (d2y/dx2)** at each point
- Uses appropriate differences based on position

### Example Output

```
Finite Differences from RK4
=========================================================
n |    x    |    y     |   y'(FD)   |  y'(exact) |   error %
=========================================================
0 |  0.0000 |  1.0000  |  0.996670  |  1.000000  |  0.3306
1 |  0.3750 |  1.2756  |  0.919521  |  0.923488  |  0.4291
2 |  0.7500 |  1.4093  |  0.574214  |  0.571476  |  0.4772
...
```

---

## Main Orchestrator (mainsolver.py)

Combines all four methods to provide comprehensive comparison.

**Workflow:**
1. Compute solutions using Euler, Heun, RK2, RK4
2. Calculate derivatives using finite differences
3. Compare against exact solution (if provided)
4. Generate comparison tables
5. Create visualization plots

**Output includes:**
- Table with y-values from all methods
- Percent relative errors vs exact solution
- Two comparison plots

### Running mainsolver.py

```bash
python mainsolver.py
```

Produces:
- Console output with formatted tables
- `Figure 1 - Numerical Methods Comparison (y-values).png` - Line plot of all methods
- `Figure 2 - Percent Relative Errors Comparison.png` - Bar chart of errors
- `output.csv` - Data export

---

## Finite Difference Module (FD.py)

Core module providing finite difference calculations.

### Compute Functions

```python
xs, ys = compute_euler()    # Euler method solution
xs, ys = compute_heun()     # Heun method solution
xs, ys = compute_rk22()     # RK2 method solution
xs, ys = compute_rk4()      # RK4 method solution
```

Each returns:
- `xs`: List of x-values
- `ys`: List of corresponding y-values

### Finite Difference Calculations

The module applies finite differences to estimate derivatives:

```python
# Forward difference at x_0
dy_dx[0] = (ys[1] - ys[0]) / h

# Central difference at interior points
for i in range(1, len(ys)-1):
    dy_dx[i] = (ys[i+1] - ys[i-1]) / (2*h)

# Backward difference at x_n
dy_dx[-1] = (ys[-1] - ys[-2]) / h

# Second derivative
d2y_dx2[i] = (ys[i+1] - 2*ys[i] + ys[i-1]) / h**2
```

---

## Error Analysis

The framework computes error metrics comparing numerical results to exact solutions.

### Percent Relative Error

```
et = |y_exact - y_numerical| / |y_exact| × 100%
```

Allows quantitative comparison of method accuracy.

### Method Comparison

- **Euler**: O(h) - fastest but least accurate
- **Heun**: O(h²) - good balance
- **RK2**: O(h²) - alternative second-order method
- **RK4**: O(h⁴) - most accurate for smooth functions

Typical accuracy order: RK4 > RK2 ≈ Heun > Euler

---

## Configuration (ivp.py)

Central configuration file for all computations.

### Parameters

```python
# Differential Equation
def f(x, y):
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0          # Starting x
y0 = 1          # Starting y
xn = 1.5        # Ending x

# Step Size
h = 0.3         # Step size

# Exact Solution (optional, for error calculation)
def y_actual_func(x):
    return x*math.cos(x) + 1

y_actual = generate_actual_solution(x0, xn, h, y_actual_func)
```

### Refinement Options

**m-Refinement (step doubling):**
```python
m = 0           # h = h_0 / 2^m
```

**n-Refinement (number of steps):**
```python
n = 10          # h = xn / (n-1)
```

---

## Utilities (utils.py)

Helper functions for display and visualization.

### Functions

```python
print_table(headers, rows)  # Pretty-print tabular data
```

Formats and displays data in aligned columns for easy reading.

### Visualization

Matplotlib plots generated automatically showing:
- All method solutions on same graph
- Percent relative error comparisons
- Exact solution overlay (if provided)

---

## Example Workflows

### Simple ODE Solution

1. Edit `ivp.py`:
```python
def f(x, y):
    return -2*x*y  # dy/dx = -2xy

x0 = 0
y0 = 1
xn = 1.0
h = 0.1
```

2. Run solver:
```bash
python mainsolver.py
```

3. View results in plots and `output.csv`

### Comparing Refinement Levels

Modify `h` in `ivp.py` and run multiple times:
```python
h = 0.3      # Run 1
h = 0.15     # Run 2
h = 0.075    # Run 3
```

Compare output to observe convergence as h decreases.

### Analyzing Derivative Approximation

Run `FD.py` to see finite difference derivatives:
```bash
python FD.py
```

Observe accuracy differences between forward, central, and backward differences.

---

## Output Files

### Console Output

**Table Format:**
```
Method Comparison
================================================
n |     x     |   y(Euler)   |   y(Heun)   |  ...
================================================
```

### Plots

**Figure 1:** y-values from all methods
- Line plot comparing solutions
- Exact solution shown if available
- Useful for visual accuracy assessment

**Figure 2:** Percent Relative Errors
- Bar chart comparing errors at each point
- Shows which method is most accurate
- Relative error = |y_exact - y_numerical| / |y_exact| × 100

### CSV Export (output.csv)

Formatted data export with:
- Step number (n)
- x-value
- y-values from each method
- Percent relative errors
- Can import into Excel/Sheets for further analysis

---

## Troubleshooting

**Issue: "No module named 'math'"**
- Solution: `import math` is in standard library, ensure Python 3.7+

**Issue: Plots not appearing**
- Solution: Install matplotlib: `pip install matplotlib`
- Or use script in interactive environment

**Issue: Large errors at later steps**
- Solution: Reduce step size (decrease h) for better accuracy
- Or use higher-order method (RK4 instead of Euler)

**Issue: Derivative values seem wrong**
- Solution: Check h is appropriate (not too large)
- Verify finite differences applied at correct positions
- Compare with exact derivative if known

---

## Dependencies

- Python 3.7+
- NumPy (if using advanced features)
- Matplotlib (for plotting)
- Math (standard library)

Install dependencies:
```bash
pip install numpy matplotlib
```

---

## Key Concepts

✓ **Finite Differences** - Approximate derivatives from discrete data  
✓ **Derivative Order** - O(h), O(h²), O(h⁴) accuracy  
✓ **Method Comparison** - Quantitative error analysis  
✓ **Error Estimation** - Percent relative error calculations  
✓ **Visualization** - Graphical method comparison  

---

## Academic Reference

This implementation demonstrates:
- Numerical solution of first-order ODEs
- Finite difference derivative approximation
- Computational accuracy and convergence
- Error analysis and comparison

Suitable for numerical analysis, differential equations, and computational methods courses.

---

## License

Academic use for ES 204 - Numerical Methods coursework.

---

**Unit 1: Numerical Differentiation** - Complete framework for solving IVPs and analyzing derivatives using finite differences.

---

## Default Problem Example

- Differential equation: `y' + x*y = cos(x)`  
- Initial condition: `y(0) = 1`  
- Interval: `0 ≤ x ≤ 2`  
- Step size: `h = 0.5`  
- Exact solution values for reference:

| x    | y_actual |
|------|----------|
| 0.0  | 1.0      |
| 0.5  | 1.32326  |
| 1.0  | 1.20079  |
| 1.5  | 0.74740  |
| 2.0  | 0.24241  |

---

This repository is designed so you can **easily replace the IVP** in `ivp.py` to solve different first-order ODEs and dynamically see the results in both table and plots.