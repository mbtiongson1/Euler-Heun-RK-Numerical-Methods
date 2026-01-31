# Euler, Heun, Runge-Kutta Numerical Methods

This project contains Python implementations of several numerical methods to solve **first-order initial value problems (IVPs)**:

- Euler Method
- Predictor–Corrector (PC / Heun)
- Ralston (RK2.2)
- Classical Runge-Kutta (RK4)

---

## How to Use

1. Modify `ivp.py` to solve your own equation:
    - Change the function `f(x, y)` to match your differential equation `y' = f(x, y)`.
    - Update `x0` and `y0` for your initial condition.
    - Adjust `h` (step size) and `xn` (end of the interval).
    - Optionally, update `y_actual` with exact values to compute percent relative error.

2. Run `mainsolver.py`:
    - This will solve the IVP using all four methods.
    - Outputs a table of y-values and percent relative errors.
    - Generates plots:
        - Figure 1: y-values of all methods compared to exact solution
        - Figure 2: Percent relative errors as a grouped bar chart
    - Automatically writes `output.csv` with all results.

3. Optional:
    - Modify `h` in `ivp.py` to refine the step size.
    - Replace `y_actual` with your exact solution if known to compute true errors.

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
