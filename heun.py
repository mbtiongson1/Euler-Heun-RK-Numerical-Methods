from ivp import f, x0, y0, h, xn#, y_actual
from utils import print_table

x = x0
yc = y0

rows = []
rows.append((0, x, "-", yc))

n = 0
while x <= xn + 1e-9:# and n < len(y_actual) - 1:
    # Predictor (Euler)
    yp = yc + h * f(x, yc)

    # Advance x
    x_next = x + h

    # Corrector (Heun)
    yc_next = yc + (h / 2) * (f(x, yc) + f(x_next, yp))

    n += 1
    x = x_next
    yc = yc_next

    # True percent relative error
    #et = (y_actual[n] - yc) / y_actual[n] * 100
    #et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, yp, yc))

print("Heun's Method / Predictor–Corrector")
print_table(
    ["n", "x", "yp (predictor)", "yc (corrected)"],
    rows
)