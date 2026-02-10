#Ralston's method

from ivp import f, x0, y0, h, xn#, y_actual
from utils import print_table

x = x0
y = y0
a1 = 1/3
a2 = 2/3
#k1 = h*f(x, y) or h*y'
#k2 = h*f(x + (3/4)*h, y+ (3/4)*k1)


rows =[]
rows.append((0, x, 0, 0, y))

num_steps = int(round((xn - x0) / h))
for n in range(1, num_steps + 1):
    k1 = h*f(x, y)
    k2 = h*f(x + (3/4)*h, y+ (3/4)*k1)
    y = y + a1*k1 + a2*k2
    x = x0 + n * h
    # True percent relative error
    #et = abs((y_actual[n] - y) / y_actual[n] * 100)
    #et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, k1, k2, y))

print("Ralston’s Method (RK2, optimal second-order)")
print_table(["n", "x", "k1", "k2", "y"], rows)