import math


# Second-order IVP form:
# y'' = f2(x, y, dy)
# where dy = y'
def f2(x, y, dy):
    return -1*dy + x - math.sin(x)


# Equivalent first-order system helpers:
# y'  = f1(x, y, dy)
# dy' = f2(x, y, dy)

# def f1(x, y, dy):
#     return dy


# Initial conditions
x0 = 0.0
y0 = 1
dy0 = -1
xn = math.pi  # limit of x. By default, starts at 0.


# Step Size, Number of steps, refinement parameters
# m = 1  # change this for m-refinement
# h = 0.3 / (2 ** m)  # m-refinement
# h = 0.3

n = 6  # change this for n-refinement
h = (xn / (n - 1))  # n-refinement


# Numerical method and polynomial degree for approximation
method = 'rk4'  # Choose based on your second-order solver (e.g., heun, rk4)
p = 10  # Polynomial degree for Vandermonde and Lagrange fits
ls_methods = ['heun', 'rk4']  # Methods to compare in least-squares fitting


# Optional exact solutions for error checking
def y_actual_func(x):
    return math.sin(x)


def dy_actual_func(x):
    return math.cos(x)


def generate_actual_solutions(x0, xn, h, y_func, dy_func):
    x = x0
    y_vals = []
    dy_vals = []
    while x <= xn + 1e-9:  # small tolerance for floating point
        y_vals.append(y_func(x))
        dy_vals.append(dy_func(x))
        x += h
    return y_vals, dy_vals


# y_actual, dy_actual = generate_actual_solutions(x0, xn, h, y_actual_func, dy_actual_func)
y_actual, dy_actual = None, None  # manual override of functions
# y_actual = [ ... ]
# dy_actual = [ ... ]
