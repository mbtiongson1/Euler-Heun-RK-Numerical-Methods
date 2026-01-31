#Ralston's method

from ivp import f, x0, y0, h, xn, y_actual
from utils import print_table

x = x0
y = y0
a1 = 1/3
a2 = 2/3
#k1 = h*f(x, y) or h*y'
#k2 = h*f(x + (3/4)*h, y+ (3/4)*k1)


rows =[]
rows.append((0, x, 0, 0, y, 0))

n = 0
while x < xn and n < len(y_actual) - 1:
    k1 = h*f(x, y)
    k2 = h*f(x + (3/4)*h, y+ (3/4)*k1)
    y = y + a1*k1 + a2*k2
    x = x + h
    n += 1
    # True percent relative error
    et = abs((y_actual[n] - y) / y_actual[n] * 100)
    et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, k1, k2, y, et_str))

print("Ralston’s Method (RK2, optimal second-order)")
print_table(["n", "x", "k1", "k2", "y", "e_t (%)"], rows)