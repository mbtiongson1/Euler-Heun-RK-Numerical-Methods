from ivp import f, x0, y0, h, xn, y_actual
from utils import print_table

x = x0
y = y0

rows =[]
rows.append((0, x, y, 0))

n = 0
while x < xn and n < len(y_actual) - 1:
    y = y + h * f(x, y)
    x = x + h
    n += 1
    # True percent relative error
    et = abs((y_actual[n] - y) / y_actual[n] * 100)
    et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, y, et_str))

print("Euler's Method")
print_table(["n", "x", "y", "e_t (%)"], rows)