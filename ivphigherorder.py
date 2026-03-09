import math

#ivp coefficients
"""
Equation:   y'' - y' + y = 0
"""
a0 = 0      #constant
a1 = 1      # y
a2 = -1     # y'
a3 = 1      # y''
a4 = 0      # y'''


# First-order
def f1(x,y):
    dy = None # = y'
    return dy

# Second-order
def f2(x, y, dy):
    d2y = -1*dy + x - math.sin(x) # = y''
    return d2y

# Third-order
def f3(x, y, dy, d2y):
    d3y = None # = y'''
    return d3y 

# Initial conditions
x0, xn = 0.0, 2.0
y0 = 0
dy0 = 3
d2y0 = 0

# Refinement
m = 1  # change this for m-refinement
h = 0.3 / (2 ** m)  # m-refinement
h = 0.3

n = 6  # change this for n-refinement
h = (xn / (n - 1))  # n-refinement

# approximating
method = 'fdm' #fdm, reduction

# fitting
p = 10  #vandermonde, lagrange
ls_methods = ['heun', 'rk4'] #least squares


# Optional exact solutions for error checking
def y_actual_func(x):
    y = None
    return y

def dy_actual_func(x):
    dy = None
    return dy

def d2y_actual_func(x):
    d2y = None
    return d2y

def generate_actual_solutions(x0, xn, h, y_func, dy_func, d2y_func):
    x = x0
    y_vals = []
    dy_vals = []
    d2y_vals = []
    while x <= xn + 1e-9:  # small tolerance for floating point
        y_vals.append(y_func(x))
        dy_vals.append(dy_func(x))
        d2y_vals.append(d2y_func(x))
        x += h
    return y_vals, dy_vals

y_actual, dy_actual, d2y_actual = generate_actual_solutions(x0, xn, h, y_actual_func, dy_actual_func, d2y_actual_func)
# y_actual, dy_actual, d2y_actual = None, None, None  # manual override of functions
# y_actual = [ ... ]
# dy_actual = [ ... ]
# d2y_actual = [ ... ]
