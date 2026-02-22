import matplotlib.pyplot as plt
import csv
from ivp import f, x0, y0, h, xn, y_actual
from utils import print_table, print_table_csv

# ---------------------------
# x_actual: independent x-coordinates for y_actual
# spaced evenly from x0 to xn based on len(y_actual)
# ---------------------------
n_actual = len(y_actual)
actual_spacing = (xn - x0) / (n_actual - 1)
x_actual = [x0 + i * actual_spacing for i in range(n_actual)]

# How many numerical steps fall between each y_actual checkpoint
actual_stride = round(actual_spacing / h)
num_steps = int(round((xn - x0) / h))

# ---------------------------
# Containers (all numerical method steps)
# ---------------------------
x_vals = [x0]

euler_y = [y0]
heun_y  = [y0]
ral_y  = [y0]
rk3_y  = [y0]
rk4_y  = [y0]

# ---------------------------
# Main loop: runs for a fixed number of steps up to xn
# ---------------------------
y_eu = y_he = y_ra = y_rk3 = y_rk = y0

for step in range(num_steps):
    x = x0 + step * h
    x_next = x0 + (step + 1) * h

    # Euler
    y_eu = y_eu + h * f(x, y_eu)

    # Heun
    yp = y_he + h * f(x, y_he)
    y_he = y_he + (h / 2) * (f(x, y_he) + f(x_next, yp))

    # Ralston
    k1 = h * f(x, y_ra)
    k2 = h * f(x + 3*h/4, y_ra + 3*k1/4)
    y_ra = y_ra + (1/3)*k1 + (2/3)*k2

    # RK3 (Kutta's 3rd order)
    k1 = h * f(x, y_rk3)
    k2 = h * f(x + h/2, y_rk3 + k1/2)
    k3 = h * f(x + h, y_rk3 - k1 + 2*k2)
    y_rk3 = y_rk3 + (1/6)*(k1 + 4*k2 + k3)

    # RK4
    k1 = h * f(x, y_rk)
    k2 = h * f(x + h/2, y_rk + k1/2)
    k3 = h * f(x + h/2, y_rk + k2/2)
    k4 = h * f(x + h,   y_rk + k3)
    y_rk = y_rk + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    x_vals.append(x_next)
    euler_y.append(y_eu)
    heun_y.append(y_he)
    ral_y.append(y_ra)
    rk3_y.append(y_rk3)
    rk4_y.append(y_rk)

# ---------------------------
# Extract numerical values at y_actual checkpoints only
# (every actual_stride steps)
# ---------------------------
euler_at_actual = [euler_y[i * actual_stride] for i in range(n_actual)]
heun_at_actual  = [heun_y [i * actual_stride] for i in range(n_actual)]
ral_at_actual   = [ral_y  [i * actual_stride] for i in range(n_actual)]
rk3_at_actual   = [rk3_y  [i * actual_stride] for i in range(n_actual)]
rk4_at_actual   = [rk4_y  [i * actual_stride] for i in range(n_actual)]

# ---------------------------
# Compute percent errors at checkpoints only
# ---------------------------
def pct_err(actual, approx):
    if actual == 0:
        return 0.0
    return abs((actual - approx) / actual * 100)

euler_e = [pct_err(y_actual[i], euler_at_actual[i]) for i in range(n_actual)]
heun_e  = [pct_err(y_actual[i], heun_at_actual[i])  for i in range(n_actual)]
ral_e   = [pct_err(y_actual[i], ral_at_actual[i])   for i in range(n_actual)]
rk3_e   = [pct_err(y_actual[i], rk3_at_actual[i])   for i in range(n_actual)]
rk4_e   = [pct_err(y_actual[i], rk4_at_actual[i])   for i in range(n_actual)]

# ---------------------------
# Prepare table rows (checkpoints only)
# ---------------------------
rows = []
for i in range(n_actual):
    rows.append((
        i, x_actual[i], y_actual[i],
        euler_at_actual[i], f"{euler_e[i]:.5f}%",
        heun_at_actual[i],  f"{heun_e[i]:.5f}%",
        ral_at_actual[i],   f"{ral_e[i]:.5f}%",
        rk3_at_actual[i],   f"{rk3_e[i]:.5f}%",
        rk4_at_actual[i],   f"{rk4_e[i]:.5f}%"
    ))

# ---------------------------
# Table headers
# ---------------------------
headers = [
    "n", "x", "y (actual)",
    "Euler (y)","Euler (error)",
    "Heun/PC (y)","Heun/PC (error)",
    "RK2.2 (y)","RK2.2 (error)",
    "RK3 (y)","RK3 (error)",
    "RK4 (y)","RK4 (error)"
]

# Print console table
print("Comparison of Numerical Methods")
print_table(headers, rows)
print_table_csv(headers, rows)

# Write CSV file
with open("output.csv", mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    for row in rows:
        formatted_row = []
        for val in row:
            if isinstance(val, str) and val.endswith('%'):
                formatted_row.append(val)
            elif isinstance(val, float):
                formatted_row.append(f"{val:.6f}")
            else:
                formatted_row.append(str(val))
        writer.writerow(formatted_row)
print("CSV file 'output.csv' successfully created!")

# ---------------------------
# Figure 1: y-values (all numerical steps + y_actual separately)
# ---------------------------
plt.figure(figsize=(10,6))
methods = [
    ("Euler",   euler_y, 'o-', 'blue'),
    ("PC", heun_y,  's-', 'green'),
    # ("RK2.2", ral_y,  '^-', 'orange'),
    # ("RK3",   rk3_y,  'p-', 'purple'),
    # ("RK4",   rk4_y,  'd-', 'red')
]

for name, y_vals, style, color in methods:
    plt.plot(x_vals, y_vals, style, label=name, color=color, linewidth=0.8,
             markersize=3)

# y_actual plotted independently at its own x-coordinates
plt.plot(x_actual, y_actual, 'k--o', label="Actual", linewidth=0.8,
         markersize=4, zorder=10)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Numerical Methods Comparison (y-values)")
plt.legend()
plt.grid(True)
plt.show()

# ---------------------------
# Figure 2: Percent error scatter (at y_actual checkpoints only)
# ---------------------------
plt.figure(figsize=(10,6))
error_methods = [
    ("Euler",   euler_e, 'o', 'blue'),
    ("PC", heun_e,  's', 'green'),
    #("RK2.2", ral_e,  '^', 'orange'),
    # ("RK3",   rk3_e,  'p', 'purple'),
    # ("RK4",   rk4_e,  'd', 'red')
]

for name, e_vals, marker, color in error_methods:
    plt.scatter(x_actual, e_vals, marker=marker, label=name, color=color, zorder=5)

plt.xlabel("x")
plt.ylabel("Percent Relative Error (%)")
plt.title("Percent Relative Errors Comparison")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
