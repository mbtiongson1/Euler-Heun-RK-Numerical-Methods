from ivpsystems import f, g, t0, x0, y0, h, tn, x_actual, y_actual
from utils import print_table

# Calculate number of steps
num_steps = int(round((tn - t0) / h))
has_actual = (
    isinstance(x_actual, list)
    and isinstance(y_actual, list)
    and len(x_actual) == num_steps + 1
    and len(y_actual) == num_steps + 1
)


def pct_err(actual, approx):
    if actual == 0:
        return 0.0
    return abs((actual - approx) / actual * 100)

rows = []
error_rows = []
if has_actual:
    ex = pct_err(x_actual[0], x0)
    ey = pct_err(y_actual[0], y0)
    rows.append((0, t0, 0, 0, 0, 0, 0, 0, 0, 0, x0, y0))
    error_rows.append((0, t0, x0, y0, x_actual[0], y_actual[0], f"{ex:.5f}%", f"{ey:.5f}%"))
else:
    rows.append((0, t0, 0, 0, 0, 0, 0, 0, 0, 0, x0, y0))

x = x0
y = y0
for n in range(1, num_steps + 1):
    t = t0 + n * h
    t_prev = t0 + (n - 1) * h

    k1x = h * f(x, y, t_prev)
    k1y = h * g(x, y, t_prev)

    k2x = h * f(x + 0.5 * k1x, y + 0.5 * k1y, t_prev + 0.5 * h)
    k2y = h * g(x + 0.5 * k1x, y + 0.5 * k1y, t_prev + 0.5 * h)

    k3x = h * f(x + 0.5 * k2x, y + 0.5 * k2y, t_prev + 0.5 * h)
    k3y = h * g(x + 0.5 * k2x, y + 0.5 * k2y, t_prev + 0.5 * h)

    k4x = h * f(x + k3x, y + k3y, t)
    k4y = h * g(x + k3x, y + k3y, t)

    x = x + (1 / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    y = y + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)

    if has_actual:
        ex = pct_err(x_actual[n], x)
        ey = pct_err(y_actual[n], y)
        rows.append((n, t, k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, x, y))
        error_rows.append((n, t, x_actual[n], y_actual[n], f"{ex:.5f}%", f"{ey:.5f}%"))
    else:
        rows.append((n, t, k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, x, y))

print("Runge-Kutta 4th Order Method (RK4) for System of IVPs")
print_table(["n", "t", "k1x", "k1y", "k2x", "k2y", "k3x", "k3y", "k4x", "k4y", "x", "y"], rows)

if has_actual:
    print("\nRK4 Accuracy Table")
    print_table(["n", "t", "x (actual)", "y (actual)", "x error", "y error"], error_rows)
