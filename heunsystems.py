from ivpsystems import f, g, t0, x0, y0, h, tn, x_actual, y_actual
from utils import print_table

t = t0
xc = x0
yc = y0

rows = []
error_rows = []
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

# First row (n=0)
if has_actual:
    ex = pct_err(x_actual[0], xc)
    ey = pct_err(y_actual[0], yc)
    rows.append((0, t, "-", "-", xc, yc))
    error_rows.append((0, t, xc, yc, x_actual[0], y_actual[0], f"{ex:.5f}%", f"{ey:.5f}%"))
else:
    rows.append((0, t, "-", "-", xc, yc))

for n in range(1, num_steps + 1):
    # Predictor (Euler)
    xp = xc + h * f(xc, yc, t)
    yp = yc + h * g(xc, yc, t)

    # Advance t
    t_next = t0 + n * h

    # Corrector PC
    xc_next = xc + (h / 2) * (f(xc, yc, t) + f(xp, yp, t_next)) #xp is xc+hf
    yc_next = yc + (h / 2) * (g(xc, yc, t) + g(xp, yp, t_next)) #xp is xc+hf

    t = t_next
    xc = xc_next
    yc = yc_next

    if has_actual:
        ex = pct_err(x_actual[n], xc)
        ey = pct_err(y_actual[n], yc)
        rows.append((n, t, xp, yp, xc, yc))
        error_rows.append((n, t, x_actual[n], y_actual[n], f"{ex:.5f}%", f"{ey:.5f}%"))
    else:
        rows.append((n, t, xp, yp, xc, yc))

print("Heun's Method for System of IVPs")
print_table(
    ["n", "t", "xp (predictor)", "yp (predictor)", "xc (corrected)", "yc (corrected)"],
    rows
)

if has_actual:
    print("\nHeun Accuracy Table")
    print_table(
        ["n", "t", "x (approx)", "y (actual)", "x error", "y error"],
        error_rows
    )
