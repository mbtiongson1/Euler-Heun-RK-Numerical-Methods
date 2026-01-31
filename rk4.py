#RK4 method

from ivp import f, x0, y0, h, xn, y_actual
from utils import print_table

x = x0
y = y0
#k1 = h*f(x, y) or h*y'
#k2 = h*f(x + (1/2)*h, y+ (1/2)*k1)
#k3 = h*f(x + (1/2)*h, y+ (1/2)*k2)
#k4 = h*f(x + h, y + k3)


rows =[]
rows.append((0, x, 0, 0, 0, 0, y, 0))

n = 0
while x < xn and n < len(y_actual) - 1:
    k1 = h*f(x, y)
    k2 = h*f(x + (1/2)*h, y+ (1/2)*k1)
    k3 = h*f(x + (1/2)*h, y+ (1/2)*k2)
    k4 = h*f(x + h, y + k3)
    y = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    x = x + h
    n += 1
    # True percent relative error
    et = abs((y_actual[n] - y) / y_actual[n] * 100)
    et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, k1, k2, k3, k4, y, et_str))

print("Runge-Kutta 4th Order Method (RK4)")
print_table(["n", "x", "k1", "k2", "k3", "k4", "y", "e_t (%)"], rows)