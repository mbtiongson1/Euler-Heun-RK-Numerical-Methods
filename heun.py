from ivp import f, x0, y0, h, xn#, y_actual
from utils import print_table

x = x0
yc = y0

rows = []
area_total = 0

# First row (n=0)
w_t = 1
area = w_t * yc
area_total += area
rows.append((0, x, "-", yc, w_t, area))
#this includes a numerical method (trapezoidal rule) to calculate the area under the curve, which is an approximation of the integral of the function. The weight w_t is used to determine how much each point contributes to the total area, with the last point having a weight of 1 and all other points having a weight of 2.
num_steps = int(round((xn - x0) / h))
for n in range(1, num_steps + 1):
    # Predictor (Euler)
    yp = yc + h * f(x, yc)

    # Advance x
    x_next = x0 + n * h

    # Corrector (Heun)
    yc_next = yc + (h / 2) * (f(x, yc) + f(x_next, yp))

    x = x_next
    yc = yc_next

    # Determine weight: 1 if last step, 2 otherwise
    w_t = 1 if n == num_steps else 2
    
    # Calculate sum product for this row
    area = w_t * yc
    area_total += area

    # True percent relative error
    #et = (y_actual[n] - yc) / y_actual[n] * 100
    #et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    rows.append((n, x, yp, yc, w_t, area))

print("Heun's Method / Predictor–Corrector")
print_table(
    ["n", "x", "yp (predictor)", "yc (corrected)", "weight", "area"],
    rows
)

# Final area total is the sum of all area values
print(f"\nArea Total: {area_total}")