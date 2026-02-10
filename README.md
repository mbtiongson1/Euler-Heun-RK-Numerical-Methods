# M-Refinement Strategy for Numerical Methods

Complete implementation of m-refinement (mesh halving) convergence study for numerical integration methods. Demonstrates how accuracy improves as step size is halved repeatedly.

## Overview

This branch focuses on **m-refinement** - a convergence strategy that reduces step size by halving (2^m), allowing analysis of method convergence rates and accuracy trends.

**Implemented Methods:**
- Euler Method
- Predictor–Corrector (Heun)
- Ralston (RK2.2)
- Classical Runge-Kutta (RK4)

**Key Feature:** m-refinement parameter in `ivp.py` allows automatic step size reduction to observe convergence behavior.

---

## Project Structure

```
├── ivp.py              # Configuration with m-refinement parameter
├── euler.py            # Euler method
├── heun.py             # Heun method
├── rk22.py             # RK2 (Ralston) method
├── rk4.py              # RK4 method
├── mainsolver.py       # Main orchestrator
├── utils.py            # Utility functions
├── test.py             # Testing framework
├── output.csv          # CSV output (generated)
└── README.md           # This file
```

---

## Quick Start

### 1. Configure Problem with m-Refinement (ivp.py)

```python
def f(x, y):
    """Your differential equation: y' = f(x, y)"""
    return math.cos(x) - x*math.sin(y)

# Initial Conditions
x0 = 0          # Starting x
y0 = 1          # Starting y
xn = 1.5        # Ending x

# M-Refinement Strategy
m = 0           # Refinement level
h = 0.3 / (2**m)  # Step size = h_0 / 2^m

# Exact Solution
def y_actual_func(x):
    return x*math.cos(x) + 1
```

### 2. Run Convergence Study

```bash
# Run with m=0
python ivp.py  # Set m = 0

# Run with m=1 (half step size)
python ivp.py  # Set m = 1

# Run with m=2 (quarter step size)
python ivp.py  # Set m = 2

# Compare results to observe convergence
python mainsolver.py
```

### 3. Observe Convergence

- **m=0:** h = 0.3 (coarsest)
- **m=1:** h = 0.15 (medium)
- **m=2:** h = 0.075 (fine)
- **m=3:** h = 0.0375 (very fine)

Each refinement level halves the step size, doubling the number of steps.

---

## M-Refinement Concept

M-refinement is a convergence study technique where step size is systematically reduced.

### Formula

```
h_m = h_0 / 2^m
```

Where:
- `h_0` = base step size
- `m` = refinement level (0, 1, 2, 3, ...)
- `h_m` = step size at refinement level m

### Example

Base step size `h_0 = 0.3`:
- m=0: h = 0.3 / 2^0 = 0.3 (10 intervals from 0 to 1.5)
- m=1: h = 0.3 / 2^1 = 0.15 (20 intervals)
- m=2: h = 0.3 / 2^2 = 0.075 (40 intervals)
- m=3: h = 0.3 / 2^3 = 0.0375 (80 intervals)

### Convergence Observation

As m increases (h decreases), numerical solutions converge to exact solution. Error typically follows:

```
Error ≈ C * h^p
```

Where `p` is the method order:
- Euler: p=1 (O(h))
- Heun/RK2: p=2 (O(h²))
- RK4: p=4 (O(h⁴))

---

## Numerical Integration Methods

### Euler Method
```
y_{n+1} = y_n + h*f(x_n, y_n)
```
- Order: O(h)
- Simplest, least accurate
- Error halves every refinement level

### Heun Method
```
Predictor:  y_p = y_n + h*f(x_n, y_n)
Corrector:  y_{n+1} = y_n + (h/2)*[f(x_n, y_n) + f(x_{n+1}, y_p)]
```
- Order: O(h²)
- Good balance
- Error quarters every refinement level

### RK2 (Ralston)
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + (2/3)*h, y_n + (2/3)*k1)
y_{n+1} = y_n + (1/4)*k1 + (3/4)*k2
```
- Order: O(h²)
- Alternative to Heun
- Similar convergence rate

### RK4 (Classical)
```
k1 = h*f(x_n, y_n)
k2 = h*f(x_n + h/2, y_n + k1/2)
k3 = h*f(x_n + h/2, y_n + k2/2)
k4 = h*f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
```
- Order: O(h⁴)
- Highly accurate
- Error reduces by factor of 16 each refinement

---

## Configuration (ivp.py)

### M-Refinement Parameter

```python
m = 0           # Change this for refinement study
h = 0.3 / (2**m)
```

### Running Refinement Study

```bash
# Manually run with different m values
for m in 0 1 2 3; do
    # Edit ivp.py, set m = $m
    python mainsolver.py
    cp output.csv output_m${m}.csv
done

# Compare output_m0.csv, output_m1.csv, etc.
```

---

## Error Convergence Analysis

### Percent Relative Error

```
e_t = |y_exact - y_numerical| / |y_exact| × 100%
```

### Convergence Rate

Expected error reduction per refinement:
- **Euler (O(h)):** Error reduces by ~2x
- **Heun/RK2 (O(h²)):** Error reduces by ~4x
- **RK4 (O(h⁴)):** Error reduces by ~16x

### Example Convergence

For RK4 solving y' = -2xy with y(0)=1 to y(1.0)=0.3679:

| m | h | Error(RK4) | Ratio |
|---|---|-----------|-------|
| 0 | 0.1 | 0.000234 | - |
| 1 | 0.05 | 0.0000146 | 16.0 |
| 2 | 0.025 | 0.00000091 | 16.0 |
| 3 | 0.0125 | 0.000000057 | 16.0 |

Notice error reduces by ~16x each level (RK4 is 4th order).

---

## Main Orchestrator (mainsolver.py)

Compares all four methods at current m-refinement level.

**Workflow:**
1. Read m from ivp.py
2. Calculate h = 0.3 / (2^m)
3. Solve using all four methods
4. Compare solutions and errors
5. Generate plots and CSV

```bash
python mainsolver.py
```

---

## Individual Method Modules

Run individual methods:
```bash
python euler.py
python heun.py
python rk22.py
python rk4.py
```

Each outputs table of solutions at current m-refinement.

---

## Example Workflows

### Convergence Study

1. Create convergence script:
```bash
for m in 0 1 2 3; do
    sed -i "s/^m = .*/m = $m/" ivp.py
    python mainsolver.py
    mv output.csv output_m${m}.csv
done
```

2. Compare results:
```bash
grep "^[0-9]*,0\." output_m*.csv | cut -d, -f1-5
```

Shows how y-values change with refinement.

### Observe Convergence Rate

Compare errors at same x-values across refinement levels to verify expected convergence rates.

### High Accuracy Solution

Set m=3 or m=4 for very accurate solution suitable for other applications.

---

## Output Files

### Console Output

```
Runge-Kutta 4th Order Method (RK4) [m=1, h=0.15]
n |        x |         y |      k1 |      k2 |      k3 |      k4
----------|---------|-----------|---------|---------|---------
0 |  0.0000  |  1.0000  |    0    |    0    |    0    |    0
...
```

Header shows m value and corresponding h.

### Plots

**Figure 1:** Solutions from all methods at current refinement level

**Figure 2:** Percent relative errors comparison

### CSV Export (output.csv)

```csv
n,x,y(Euler),y(Heun),y(RK2),y(RK4),y_exact,error_Euler,...
0,0.0000,1.0000,1.0000,1.0000,1.0000,1.0000,0.0000,...
```

Includes m value in output for tracking.

---

## Utilities (utils.py)

```python
print_table(headers, rows)  # Pretty-print data
```

---

## Troubleshooting

**Issue: Results don't change between refinement levels**
- Solution: Verify m value changed in ivp.py
- Check h = 0.3 / (2**m) is recalculated

**Issue: Method becomes unstable at high m**
- Solution: Some problems require special handling
- Check if problem is stiff

**Issue: Convergence rate doesn't match theory**
- Solution: May indicate implementation issue
- Or exact solution function is wrong
- Check problem parameters

---

## Dependencies

- Python 3.7+
- Matplotlib
- Math (standard library)

---

## Key Concepts

✓ **M-Refinement** - Step doubling convergence study  
✓ **Convergence Rate** - O(h^p) behavior  
✓ **Error Reduction** - Systematic error analysis  
✓ **Method Comparison** - Accuracy across refinement levels  

---

## Academic Reference

Demonstrates:
- Convergence analysis of numerical methods
- Error estimation and reduction
- Method order verification
- Practical convergence rate calculation

---

## License

Academic use for ES 204 - Numerical Methods coursework.

---

**Unit 1: M-Refinement Strategy** - Convergence study of numerical integration methods through systematic step size reduction.

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