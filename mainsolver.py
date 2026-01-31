import matplotlib.pyplot as plt
import csv
import numpy as np
from ivp import f, x0, y0, h, xn, y_actual
from utils import print_table, print_table_csv

# ---------------------------
# Containers
# ---------------------------
x_vals = [x0]

euler_y, euler_e = [y0], [0.0]
heun_y,  heun_e  = [y0], [0.0]
ral_y,   ral_e   = [y0], [0.0]
rk4_y,   rk4_e   = [y0], [0.0]

# ---------------------------
# Main loop for all methods
# ---------------------------
x = x0
y_eu = y_he = y_ra = y_rk = y0
n = 0

while x < xn and n < len(y_actual) - 1:
    x_next = x + h

    # Euler
    y_eu = y_eu + h * f(x, y_eu)

    # Heun
    yp = y_he + h * f(x, y_he)
    y_he = y_he + (h / 2) * (f(x, y_he) + f(x_next, yp))

    # Ralston
    k1 = h * f(x, y_ra)
    k2 = h * f(x + 3*h/4, y_ra + 3*k1/4)
    y_ra = y_ra + (1/3)*k1 + (2/3)*k2

    # RK4
    k1 = h * f(x, y_rk)
    k2 = h * f(x + h/2, y_rk + k1/2)
    k3 = h * f(x + h/2, y_rk + k2/2)
    k4 = h * f(x + h,   y_rk + k3)
    y_rk = y_rk + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    # Advance
    x = x_next
    n += 1
    x_vals.append(x)

    # Store values
    euler_y.append(y_eu)
    heun_y.append(y_he)
    ral_y.append(y_ra)
    rk4_y.append(y_rk)

    # Compute true percent relative errors
    euler_e.append(abs((y_actual[n] - y_eu) / y_actual[n] * 100))
    heun_e.append(abs((y_actual[n] - y_he) / y_actual[n] * 100))
    ral_e.append(abs((y_actual[n] - y_ra) / y_actual[n] * 100))
    rk4_e.append(abs((y_actual[n] - y_rk) / y_actual[n] * 100))

# ---------------------------
# Prepare table rows
# ---------------------------
rows = []
for i in range(len(x_vals)):
    rows.append((
        i, x_vals[i], y_actual[i],
        euler_y[i], f"{euler_e[i]:.5f}%",
        heun_y[i],  f"{heun_e[i]:.5f}%",
        ral_y[i],   f"{ral_e[i]:.5f}%",
        rk4_y[i],   f"{rk4_e[i]:.5f}%"
    ))

# ---------------------------
# Table headers
# ---------------------------
headers = [
    "n", "x", "y (actual)",
    "Euler (y)","Euler (error)",
    "Heun/PC (y)","Heun/PC (error)",
    "RK2.2 (y)","RK2.2 (error)",
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
# Figure 1: y-values with non-overlapping labels
# ---------------------------
plt.figure(figsize=(10,6))
methods = [
    ("Euler", euler_y, 'o-', 'blue'),
    ("Heun/PC", heun_y,  's-', 'green'),
    ("RK2.2", ral_y, '^-', 'orange'),
    ("RK4", rk4_y,   'd-', 'red')
]

for idx, (name, y_vals, style, color) in enumerate(methods):
    plt.plot(x_vals, y_vals, style, label=name, color=color)
    
    # Label the first point only once
    if idx == 0:
        plt.annotate(
            f"{y_vals[0]:.5f}",
            xy=(x_vals[0], y_vals[0]),
            xytext=(0, 15),  # offset above
            textcoords='offset points',
            arrowprops=dict(arrowstyle='->', color='black', lw=0.8),
            fontsize=7,
            ha='center'
        )
    
    # Label the rest of the points, skipping n=0
    for n, (xi, yi) in enumerate(zip(x_vals[1:], y_vals[1:]), start=1):
        dx = 5 * ((idx+1) * n)
        dy = 6 * (idx+1)
        plt.annotate(
            f"{yi:.5f}",
            xy=(xi, yi),
            xytext=(dx, dy),
            textcoords='offset points',
            arrowprops=dict(arrowstyle='->', color=color, lw=0.8),
            fontsize=7,
            ha='center'
        )

# Actual solution
plt.plot(x_vals, y_actual[:len(x_vals)], 'k--', label="Actual")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Numerical Methods Comparison (y-values)")
plt.legend()
plt.grid(True)
plt.subplots_adjust(bottom=0.35)  # margin for labels
plt.show()

# ---------------------------
# Figure 2: Percent error bar chart with top labels
# ---------------------------
plt.figure(figsize=(10,6))
x_indices = np.arange(len(x_vals))
width = 0.2

bars = [
    plt.bar(x_indices - 1.5*width, euler_e, width, label="Euler", color='blue'),
    plt.bar(x_indices - 0.5*width, heun_e, width, label="Heun/PC", color='green'),
    plt.bar(x_indices + 0.5*width, ral_e, width, label="RK2.2", color='orange'),
    plt.bar(x_indices + 1.5*width, rk4_e, width, label="RK4", color='red')
]

# Add data labels on top of each bar
for idx, bar_group in enumerate(bars):
    for n, bar in enumerate(bar_group):
        height = bar.get_height()
        # Only label first group once (n=0), then label the rest normally
        if n == 0 and idx == 0:
            plt.text(
                bar.get_x() + bar.get_width()/2,
                height + 0.5,
                "0%",
                ha='center',
                va='bottom',
                fontsize=8
            )
        elif n != 0:
            plt.text(
                bar.get_x() + bar.get_width()/2,
                height + 0.5,
                f"{height:.2f}%",
                ha='center',
                va='bottom',
                fontsize=8
            )

plt.xticks(x_indices, [f"{x:.2f}" for x in x_vals])
plt.xlabel("x")
plt.ylabel("Percent Relative Error (%)")
plt.title("Percent Relative Errors Comparison")
plt.legend()
plt.grid(axis='y')
plt.show()