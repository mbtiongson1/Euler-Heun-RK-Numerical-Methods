from ivphigherorder import f2, x0, y0, dy0, d2y0, h, xn, a0, a1, a2, a3, a4, y_actual
from matrix import matrixsolver
import math
import numpy as np
import pandas as pd

print(f"GE: {a3}y'' + {a2}y' + {a1}y = {a0}")

# Discretization
def discrete_domain():
    x = x0
    ystr = ""
    x_i = []
    y_i = []
    i = 0
    while x <= xn + 1e-9:
        ystr = "y"+str(i)
        y_i.append(ystr)
        x_i.append(x)
        i += 1
        x += h
    print(f'Domain discretized:')
    print(f'x_i = {x_i}')
    print(f'y_i = {y_i}')
    return x_i, y_i

# CD methods
def centraldiff1(y_before, y_after):
    y_b, y_a = y_before, y_after
    dy = (y_b - y_a) / (2*h)
    return dy

def centraldiff2(y_before, y_current, y_after):
    y_b, y_c, y_a = y_before, y_current, y_after
    d2y = (y_a - 2*y_c + y_b) / (h**2)
    return d2y

#coefficients
def discrete_ge():
    #a_after
    a_after_d2y = a3 * 2
    a_after_dy = a2 * h
    a_after = a_after_d2y + a_after_dy
    # a_current
    a_current_d2y = a3 * -2 * 2
    a_current_dy = 0
    a_current_y = a1 * 2 * h**2
    a_current = a_current_d2y + a_current_dy + a_current_y
    #a_before
    a_before_d2y = a3 * 2
    a_before_dy = a2 * -1 * h
    a_before = a_before_d2y + a_before_dy

    c = 2*(h**2) * a0
    """
    Equation: a*y_{i-1} + a*y_i + a*y_{i+1} = c
    """
    print(f'\nSystem of Equations:')
    print(f'{a_before}y_(i-1) + {a_current}y_i + {a_after}y_(i+1) = {c}')
    return a_before, a_current, a_after, c

def option(x, constant): #requires the constant term
    if x==3:
        a_current = 4
        a_after = -1
        a_before = 0
        c = constant * 2*h + (3*y0) + dy0
        """
        Equation: (-3*y0 + 4*y1 - y2) / 2*h = dy0
        """
        print(f'{a_before}y0 + {a_current}y1 + {a_after}y2 = {c}')
        return a_current, a_after, c
    else:
        return None, None, None

def errors(y_list):
    y = y_list
    if len(y) == len(y_actual):
        error_list = []
        for i, y_val in enumerate(y):
            if y_actual[i] == 0:
                error_calc = 0
            else:
                error_calc = abs(y_val - y_actual[i]) / (y_actual[i])
            error_list.append(error_calc)
        return error_list
    else: 
        return None


def printtable(x_i, y_i, y_actual, error_i):
    if y_actual is None:
        y_actual = [None] * len(y_i)
    if error_i is None:
        error_i = [None] * len(y_i)

    table = pd.DataFrame(
        {
            "x_i": x_i,
            "y_i": y_i,
            "y_actual": y_actual,
            "error_i": error_i,
        }
    )
    print(table.to_string(index=False))


x_i = []
y_i = []
x_i, y_i = discrete_domain()
y_i = [y0] #reinitializing after discretizing

#getting constants
a_b, a_c, a_a, c = discrete_ge()
# a_c1, a_a1, c1 = option(1, c)
# a_c2, a_a2, c2 = option(2, c)
# print(f"constant: {c}")
a_c3, a_a3, c3 = option(3, a0)

#for matrix
A = []
row1 = []
row2 = []
row3 = []
row4 = []       #for 4x4
# row5 = []     #for 5x5

#A y1, y2, y3, y4 = y0
row1 = [a_c3, a_a3, 0, 0] #option3
row2 = [a_c, a_a, 0, 0]
row3 = [a_b, a_c, a_a, 0]
row4 = [0, a_b, a_c, a_a]

A = [
    row1, 
    row2, 
    row3, 
    row4
    ]
# row5 = [0, 0, a_b, a_c, a_a]

#C constants vector
C = [c3 + y0,
    a_b * y0,
    c,
    c #add one more c for 5x5
    ]

print("[A] = ")
A_np = np.array(A, dtype=float)
print(f'\n {A_np}')

print(C)

y_coeffs = matrixsolver(A, C)

y_i.extend(y_coeffs.tolist())
error_i = errors(y_i)
# print(f'y_values: {y_i}')
# print(f'errors: {error_i}')
printtable(x_i, y_i, y_actual, error_i)
