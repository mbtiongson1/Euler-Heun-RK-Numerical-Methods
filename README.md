# Numerical Integration Methods

Complete implementation of numerical methods for solving first-order initial value problems (IVPs) using Euler, Heun, Runge-Kutta, and related integration schemes.

## Overview

This branch focuses on **numerical integration** - solving differential equations using step-by-step approximation methods.

**Implemented Methods:**
- Euler Method
- Predictor–Corrector (Heun)
- Ralston (RK2.2)
- Classical Runge-Kutta (RK4)

---

## Project Structure

```
├── ivp.py              # Central configuration file
├── euler.py            # Euler method implementation
├── heun.py             # Heun method implementation
├── rk22.py             # RK2 (Ralston) implementation
├── rk4.py              # RK4 method implementation
├── mainsolver.py       # Main orchestrator
├── utils.py            # Utility functions
├── test.py             # Testing framework
├── output.csv          # CSV output (generated)
└── README.md           # This file
```

---

## Quick Start

### 1. Configure Your Problem (ivp.py)

```python
def f(x, y):
    """Your differential equation: y' = f(x, y)"""
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0          # Starting x
y0 = 1          # Starting y
xn = 1.5        # Ending x

# Step Size
h = 0.3         # Fixed step size

# Exact Solution (optional, for error calculation)
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
```

### 3. View Results

- **Console Output:** Tables with y-values and errors
- **Plots:** Comparison figures
- **CSV Export:** `output.csv` with detailed data

---

## Numerical Integration Methods

All methods solve the initial value problem:
```
dy/dx = f(x, y)
y(x0) = y0
```

### Euler Method
**First-order method** - simplest but least accurate.
```
y_{n+1} = y_n + h*f(x_n, y_n)
```
- Local error: O(h²)
- Global error: O(h)
- Good for understanding, poor for accuracy

### Heun Method (Improved Euler)
**Second-order predictor-corrector method** - two function evaluations per step.
```
Predictor:  y_p = y_n + h*f(x_n, y_n)
Corrector:  y_{n+1} = y_n + (h/2)*[f(x_n, y_n) + f(x_{n+1}, y_p)]
```
- Local error: O(h³)
- Global error: O(h²)
- Good balance of accuracy and simplicity

### RK2 (Ralston Method)
**Alternative second-order method** - uses weighted evaluations.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + (2/3)*h, y_n + (2/3)*k1)
y_{n+1} = y_n + (1/4)*k1 + (3/4)*k2
```
- Local error: O(h³)
- Global error: O(h²)
- Similar accuracy to Heun, different formula

### RK4 (Classical Runge-Kutta)
**Fourth-order method** - most widely used, excellent accuracy.
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h/2, y_n + k2/2)
k4 = h*f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
```
- Local error: O(h⁵)
- Global error: O(h⁴)
- Standard choice for smooth functions

---

## Configuration (ivp.py)

Central configuration for all computations.

### Parameters

```python
# Differential Equation
def f(x, y):
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0          # Initial x value
y0 = 1          # Initial y value
xn = 1.5        # End point

# Step Size
h = 0.3         # Step size (or use refinement)

# Exact Solution (optional)
def y_actual_func(x):
    return x*math.cos(x) + 1
```

### Refinement Options

**Fixed Step Size:**
```python
h = 0.3         # Direct specification
```

**m-Refinement (halving):**
```python
m = 0           # h = h_0 / 2^m (doubling number of steps)
h = 0.3 / (2**m)
```

**n-Refinement (number of points):**
```python
n = 10          # h = xn / (n-1)
h = xn / (n - 1)
```

---

## Main Orchestrator (mainsolver.py)

Solves the IVP using all four methods and compares results.

**Workflow:**
1. Compute solutions using Euler, Heun, RK2, RK4
2. Compare solutions at each step
3. Calculate percent relative errors vs exact solution
4. Generate comparison tables and plots
5. Export results to CSV

**Output includes:**
- Console table with y-values from all methods
- Percent relative error for each method
- Two matplotlib figures
- CSV file with all data

### Running mainsolver.py

```bash
python mainsolver.py
```

Generates:
- `Figure 1 - Numerical Methods Comparison (y-values).png`
- `Figure 2 - Percent Relative Errors Comparison.png`
- `output.csv` - Data export

---

## Individual Method Modules

Each method can be run independently:

### euler.py
```bash
python euler.py
```
Outputs Euler method results in table format.

### heun.py
```bash
python heun.py
```
Outputs Heun method results with predictor-corrector values.

### rk22.py
```bash
python rk22.py
```
Outputs RK2 (Ralston) method results.

### rk4.py
```bash
python rk4.py
```
Outputs RK4 method results.

---

## Error Analysis

Compares numerical solutions to exact solution.

### Percent Relative Error

```
e_t = |y_exact - y_numerical| / |y_exact| × 100%
```

Shows accuracy of each method.

### Typical Error Order

From most to least accurate:
1. **RK4** - O(h⁴) - Best for smooth functions
2. **RK2 / Heun** - O(h²) - Good balance
3. **Euler** - O(h) - Least accurate

---

## Utilities (utils.py)

Helper functions for output and visualization.

### print_table(headers, rows)

Pretty-prints tabular data with aligned columns.

```python
from utils import print_table

print_table(
    ["n", "x", "y(Euler)", "y(Heun)", "y(RK4)"],
    rows
)
```

---

## Example Workflows

### Basic Usage

1. Set up ODE in `ivp.py`:
```python
def f(x, y):
    return -2*x*y  # dy/dx = -2xy

x0 = 0
y0 = 1
xn = 1.0
h = 0.1
```

2. Run:
```bash
python mainsolver.py
```

3. View results in plots and CSV

### Convergence Study

Test multiple step sizes:
```python
h = 0.1      # Run 1
h = 0.05     # Run 2  
h = 0.025    # Run 3
```

Compare outputs to observe convergence behavior.

### Single Method

Solve with just one method:
```bash
python rk4.py
```

---

## Output Files

### Console Output

```
Runge-Kutta 4th Order Method (RK4)
n |        x |         y |      k1 |      k2 |      k3 |      k4
----------|---------|-----------|---------|---------|---------
0 |  0.0000  |  1.0000  |    0    |    0    |    0    |    0
1 |  0.1000  |  1.0997  | 0.0997  | 0.0947  | 0.0948  | 0.0900
...
```

### Plots

**Figure 1:** y-values comparison
- Line plot showing all method solutions
- Exact solution overlay
- Visual accuracy comparison

**Figure 2:** Percent relative errors
- Bar chart comparing errors
- Shows which method is most accurate at each point
- Error magnitude visualization

### CSV Export (output.csv)

```csv
n,x,y(Euler),y(Heun),y(RK2),y(RK4),y_exact,error_Euler,error_Heun,...
0,0.0000,1.0000,1.0000,1.0000,1.0000,1.0000,0.0000,0.0000,...
```

---

## Troubleshooting

**Issue: Large errors with small h**
- Solution: Check `y_actual_func` is correct
- Or remove error calculation if exact solution unknown

**Issue: Plots not showing**
- Solution: Install matplotlib: `pip install matplotlib`

**Issue: Method diverges**
- Solution: Check f(x,y) is correctly implemented
- Reduce h or use higher-order method (RK4)

**Issue: Overflow errors**
- Solution: Problem may be stiff or unstable
- Try much smaller h or use different method

---

## Dependencies

- Python 3.7+
- Matplotlib (for plotting)
- Math (standard library)

Install:
```bash
pip install matplotlib
```

---

## Key Concepts

✓ **Initial Value Problems** - Solving ODE with known initial condition  
✓ **Numerical Integration** - Step-by-step approximation of solution  
✓ **Error Estimation** - Local and global truncation errors  
✓ **Method Comparison** - Accuracy vs complexity trade-offs  
✓ **Convergence** - Behavior as h → 0  

---

## Academic Reference

This implementation demonstrates:
- Four major numerical integration methods
- Error analysis and convergence
- Comparative method evaluation
- Computational accuracy assessment

Suitable for numerical analysis, differential equations, and computational methods courses.

---

## License

Academic use for ES 204 - Numerical Methods coursework.

---

**Unit 1: Numerical Integration** - Solving initial value problems using Euler, Heun, RK2, and RK4 methods.