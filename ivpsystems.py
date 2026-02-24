import math


# System of IVPs with independent variable t
# x' = f(x, y, t)
# y' = g(x, y, t)
def f(x, y, t):  # x' define this!
   return -0.25*x

def g(x, y, t):  # y' define this!
    return 5 + 0.1*x - 0.4*y


# Initial conditions
t0 = 0
x0 = 1
y0 = -2
tn = 2  # limit of t. By default, starts at 0.

# Step Size, Number of steps, refinement parameters
m = 0  # change this for m-refinement
h = 0.5 / (2 ** m)  # m-refinement
# h = 0.3

# n = 9  # change this for n-refinement
# h = (tn / (n - 1))  # n-refinement


# Numerical method and polynomial degree for approximation
method = 'rk4'  # Choose: heun, rk4
p = 10  # Polynomial degree for Vandermonde and Lagrange fits
ls_methods = ['heun', 'rk4']  # Methods to compare in least-squares fitting


# Exact solutions for error checking
def x_exact_func(t):
    return 4 * math.exp(-0.5 * t) #define this if given


def y_exact_func(t):
    return 2 * math.exp(-1.1 * t) * (
        math.exp(0.6 * t) - 4.6667 * math.exp(0.8 * t) + 0.66667 * math.exp(1.1 * t) #define this if given
    )


def generate_actual_solutions(t0, tn, h, x_func, y_func):
    t = t0
    x_vals = []
    y_vals = []
    while t <= tn + 1e-9:  # small tolerance for floating point
        x_vals.append(x_func(t))
        y_vals.append(y_func(t))
        t += h
    return x_vals, y_vals


x_actual, y_actual = generate_actual_solutions(t0, tn, h, x_exact_func, y_exact_func)

#manual override of functions
x_actual = [1, 0.882497, 0.78801, 0.687289, 0.606531]
y_actual = [-2, 0.670915, 2.852668, 4.3455, 6.08953]