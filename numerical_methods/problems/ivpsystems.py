import math


# System of IVPs with independent variable x
# y' = f(y, z, x)
# z' = g(y, z, x)
def f(y, z, x):  # z=y' define this!
   return z

def g(y, z, x):  # z'=y'' define this!
    return -1*y + x**2 - 2*math.sin(x)


# Initial conditions if x,y,t
# x0 = 1
# y0 = -1
# t0 = 0
# tn = 0.2  # limit of t. By default, starts at 0.

#Initial conditions if y,z,x
x0 = 0
xn = 1
y0 = 1

beta = -1.57119 #given as y'10

def linearize(gamma1, gamma2, g1, g2, beta):
    if (g1-g2) == 0: 
        return 0
    else:
        gamma3 = gamma2 - (g2-beta) * (gamma2-gamma1)/(g2-g1)
        return gamma3

gamma1 = 3
gamma2 = -3
g1 = 0
g2 = 0

z0 = gamma1
gamma3 = linearize(gamma1, gamma2, g1, g2, beta)

# Normalizing:
t0 = x0
tn = xn
x0 = y0
y0 = z0

# Step Size, Number of steps, refinement parameters
m = 0  # change this for m-refinement
h = 0.05 / (2 ** m)  # m-refinement
h = 0.1 # override

# n = 9  # change this for n-refinement
# h = (tn / (n - 1))  # n-refinement


# Numerical method and polynomial degree for approximation
method = 'rk4'  # Choose: heun, rk4
p = 10  # Polynomial degree for Vandermonde and Lagrange fits
ls_methods = ['heun', 'rk4']  # Methods to compare in least-squares fitting


# Exact solutions for error checking
def x_exact_func(t):
    return None #define this if given


def y_exact_func(t):
    return None

def generate_actual_solutions(t0, tn, h, x_func, y_func):
    t = t0
    x_vals = []
    y_vals = []
    while t <= tn + 1e-9:  # small tolerance for floating point
        x_vals.append(x_func(t))
        y_vals.append(y_func(t))
        t += h
    return x_vals, y_vals


# x_actual, y_actual = generate_actual_solutions(t0, tn, h, x_exact_func, y_exact_func)
x_actual, y_actual = None, None # manual override of functions
x_actual = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
y_actual = [1.0, 0.95675, ]