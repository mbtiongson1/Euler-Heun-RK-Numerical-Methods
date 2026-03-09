from ivphigherorder import f2, x0, y0, dy0, d2y0, h, xn, a0, a1, a2, a3, a4
from matrix import solve_linear_system
from utils import print_table

# Discretization
def discrete_domain():
    x = x0
    ystr = ""
    x_i = [x0]
    y_i = [y0]
    i = 0
    while x <= xn + 1e-9:
        ystr = "y"+str(i)
        y_i.append(ystr)
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
    #a_before
    a_before_d2y, a_after_d2y = a3 * 2
    a_before_dy = a2 * h
    a_before = a_before_d2y + a_before_dy
    # a_current
    a_current_d2y = a3 * -2 * 2
    a_current_dy = 0
    a_current = a_current_d2y + a_current_dy
    #a_after
    a_after_dy = -1 * a2 * h
    a_after = a_after_d2y + a_after_dy

    c = 2*(h**2) * a0
    """
    Equation: a*y_{i-1} + a*y_i + a*y_{i+1} = c
    """
    return a_before, a_current, a_after, c

def option3(constant): #requires the constant term
    a_current = 4
    a_after = -1
    c = constant * 2*h + 3*y0
    """
    Equation: (-3*y0 + 4*y1 - y2) / 2*h = c
    """
    return a_current, a_after, c

x_i = []
y_i = []
x_i, y_i = discrete_domain()

#getting constants
a_b, a_c, a_a, c1 = discrete_ge()
a_c2, a_a2, c2 = option3(c1)

#for matrix
A = []
row1 = []
row2 = []
row3 = []
row4 = []       #for 4x4
# row5 = []     #for 5x5

# for row, x in enumerate(x_i):
    
