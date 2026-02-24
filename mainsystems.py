import csv

from ivpsystems import f, g, t0, x0, y0, h, tn, x_actual, y_actual
from utils import print_table, print_table_csv

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False

# Choose any subset from: euler, heun, ral, rk3, rk4
# active_methods = ['euler', 'heun', 'ral', 'rk3', 'rk4']
active_methods = ['rk4']
# ---------------------------
# Time grid
# ---------------------------
num_steps = int(round((tn - t0) / h))
t_vals = [t0]

# ---------------------------
# Containers (all numerical method steps)
# ---------------------------
euler_x = [x0]
euler_y = [y0]

heun_x = [x0]
heun_y = [y0]

ral_x = [x0]
ral_y = [y0]

rk3_x = [x0]
rk3_y = [y0]

rk4_x = [x0]
rk4_y = [y0]

# ---------------------------
# Main loop
# ---------------------------
x_eu = x_he = x_ra = x_rk3 = x_rk4 = x0
y_eu = y_he = y_ra = y_rk3 = y_rk4 = y0

for step in range(num_steps):
    t = t0 + step * h
    t_next = t0 + (step + 1) * h

    # Euler
    kx = h * f(x_eu, y_eu, t)
    ky = h * g(x_eu, y_eu, t)
    x_eu = x_eu + kx
    y_eu = y_eu + ky

    # Heun
    xp = x_he + h * f(x_he, y_he, t)
    yp = y_he + h * g(x_he, y_he, t)
    x_he = x_he + (h / 2) * (f(x_he, y_he, t) + f(xp, yp, t_next))
    y_he = y_he + (h / 2) * (g(x_he, y_he, t) + g(xp, yp, t_next))

    # Ralston (RK2.2)
    k1x = h * f(x_ra, y_ra, t)
    k1y = h * g(x_ra, y_ra, t)
    k2x = h * f(x_ra + 3 * k1x / 4, y_ra + 3 * k1y / 4, t + 3 * h / 4)
    k2y = h * g(x_ra + 3 * k1x / 4, y_ra + 3 * k1y / 4, t + 3 * h / 4)
    x_ra = x_ra + (1 / 3) * k1x + (2 / 3) * k2x
    y_ra = y_ra + (1 / 3) * k1y + (2 / 3) * k2y

    # RK3 (Kutta's 3rd order)
    k1x = h * f(x_rk3, y_rk3, t)
    k1y = h * g(x_rk3, y_rk3, t)
    k2x = h * f(x_rk3 + k1x / 2, y_rk3 + k1y / 2, t + h / 2)
    k2y = h * g(x_rk3 + k1x / 2, y_rk3 + k1y / 2, t + h / 2)
    k3x = h * f(x_rk3 - k1x + 2 * k2x, y_rk3 - k1y + 2 * k2y, t + h)
    k3y = h * g(x_rk3 - k1x + 2 * k2x, y_rk3 - k1y + 2 * k2y, t + h)
    x_rk3 = x_rk3 + (1 / 6) * (k1x + 4 * k2x + k3x)
    y_rk3 = y_rk3 + (1 / 6) * (k1y + 4 * k2y + k3y)

    # RK4
    k1x = h * f(x_rk4, y_rk4, t)
    k1y = h * g(x_rk4, y_rk4, t)
    k2x = h * f(x_rk4 + k1x / 2, y_rk4 + k1y / 2, t + h / 2)
    k2y = h * g(x_rk4 + k1x / 2, y_rk4 + k1y / 2, t + h / 2)
    k3x = h * f(x_rk4 + k2x / 2, y_rk4 + k2y / 2, t + h / 2)
    k3y = h * g(x_rk4 + k2x / 2, y_rk4 + k2y / 2, t + h / 2)
    k4x = h * f(x_rk4 + k3x, y_rk4 + k3y, t_next)
    k4y = h * g(x_rk4 + k3x, y_rk4 + k3y, t_next)
    x_rk4 = x_rk4 + (1 / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    y_rk4 = y_rk4 + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)

    t_vals.append(t_next)
    euler_x.append(x_eu)
    euler_y.append(y_eu)
    heun_x.append(x_he)
    heun_y.append(y_he)
    ral_x.append(x_ra)
    ral_y.append(y_ra)
    rk3_x.append(x_rk3)
    rk3_y.append(y_rk3)
    rk4_x.append(x_rk4)
    rk4_y.append(y_rk4)


def pct_err(actual, approx):
    if actual == 0:
        return 0.0
    return abs((actual - approx) / actual * 100)


series_x = {
    "euler": euler_x,
    "heun": heun_x,
    "ral": ral_x,
    "rk3": rk3_x,
    "rk4": rk4_x,
}

series_y = {
    "euler": euler_y,
    "heun": heun_y,
    "ral": ral_y,
    "rk3": rk3_y,
    "rk4": rk4_y,
}

method_labels = {
    "euler": "Euler",
    "heun": "PC",
    "ral": "RK2.2",
    "rk3": "RK3",
    "rk4": "RK4",
}

# ---------------------------
# Checkpoint selection (aligned to actual solution if available)
# ---------------------------
has_actual = (
    isinstance(x_actual, list)
    and isinstance(y_actual, list)
    and len(x_actual) > 1
    and len(x_actual) == len(y_actual)
)

if has_actual:
    n_actual = len(x_actual)
    actual_spacing = (tn - t0) / (n_actual - 1)
    actual_stride = max(1, round(actual_spacing / h))
    t_check = [t0 + i * actual_spacing for i in range(n_actual)]
    indices = [i * actual_stride for i in range(n_actual)]
else:
    t_check = t_vals
    indices = list(range(len(t_vals)))

# Step: validate whether actual arrays are truly usable for plotting.
show_actual_in_plots = (
    has_actual
    and all(isinstance(v, (int, float)) for v in x_actual)
    and all(isinstance(v, (int, float)) for v in y_actual)
    and len(t_check) == len(x_actual) == len(y_actual)
)


# ---------------------------
# Tables
# ---------------------------
if has_actual:
    rows_x = []
    rows_y = []
    for i, idx in enumerate(indices):
        row_x = [i, t_check[i], x_actual[i]]
        row_y = [i, t_check[i], y_actual[i]]
        for m in active_methods:
            approx_x = series_x[m][idx]
            approx_y = series_y[m][idx]
            row_x.extend([approx_x, f"{pct_err(x_actual[i], approx_x):.5f}%"])
            row_y.extend([approx_y, f"{pct_err(y_actual[i], approx_y):.5f}%"])
        rows_x.append(tuple(row_x))
        rows_y.append(tuple(row_y))

    headers_x = ["n", "t", "x (actual)"]
    headers_y = ["n", "t", "y (actual)"]
    for m in active_methods:
        label = method_labels[m]
        headers_x.extend([f"{label} (x)", f"{label} (error)"])
        headers_y.extend([f"{label} (y)", f"{label} (error)"])

    print("Comparison of Numerical Methods for x(t)")
    print_table(headers_x, rows_x)
    print_table_csv(headers_x, rows_x)

    print("\nComparison of Numerical Methods for y(t)")
    print_table(headers_y, rows_y)
    print_table_csv(headers_y, rows_y)

    with open("output_systems_x.csv", mode="w", newline="") as fx:
        writer = csv.writer(fx)
        writer.writerow(headers_x)
        writer.writerows(rows_x)

    with open("output_systems_y.csv", mode="w", newline="") as fy:
        writer = csv.writer(fy)
        writer.writerow(headers_y)
        writer.writerows(rows_y)

    print("\nCSV files 'output_systems_x.csv' and 'output_systems_y.csv' successfully created!")
else:
    headers = ["n", "t"]
    for m in active_methods:
        label = method_labels[m]
        headers.extend([f"{label} (x)", f"{label} (y)"])

    rows = []
    for i, t in enumerate(t_vals):
        row = [i, t]
        for m in active_methods:
            row.extend([series_x[m][i], series_y[m][i]])
        rows.append(tuple(row))

    print("Comparison of Numerical Methods (no actual solutions provided)")
    print_table(headers, rows)
    print_table_csv(headers, rows)


if HAS_MATPLOTLIB:
    # ---------------------------
    # Plots
    # ---------------------------
    styles = {
        "euler": ('o-', 'blue'),
        "heun": ('s-', 'green'),
        "ral": ('^-', 'orange'),
        "rk3": ('p-', 'purple'),
        "rk4": ('d-', 'red'),
    }

    x_palette = ['#1f77b4', '#2a9df4', '#1769aa', '#3b82f6', '#0ea5e9']
    y_palette = ['#d62728', '#ef4444', '#b91c1c', '#f97316', '#fb7185']
    label_stride = max(1, len(t_vals) // 12)

    plt.figure(figsize=(12, 7))
    for idx, m in enumerate(active_methods):
        x_color = x_palette[idx % len(x_palette)]
        y_color = y_palette[idx % len(y_palette)]

        plt.plot(
            t_vals, series_x[m], '-o',
            label=f"{method_labels[m]} x(t)",
            color=x_color, linewidth=1.0, markersize=3
        )
        plt.plot(
            t_vals, series_y[m], '--s',
            label=f"{method_labels[m]} y(t)",
            color=y_color, linewidth=1.0, markersize=3
        )

        for i in range(0, len(t_vals), label_stride):
            plt.annotate(
                f"{series_x[m][i]:.2f}",
                xy=(t_vals[i], series_x[m][i]),
                xytext=(0, 6),
                textcoords='offset points',
                fontsize=7,
                color=x_color,
                ha='center'
            )
            plt.annotate(
                f"{series_y[m][i]:.2f}",
                xy=(t_vals[i], series_y[m][i]),
                xytext=(0, -10),
                textcoords='offset points',
                fontsize=7,
                color=y_color,
                ha='center'
            )

    if show_actual_in_plots:
        plt.plot(t_check, x_actual, 'k-o', label="Actual x(t)", linewidth=1.2, markersize=3, zorder=10)
        plt.plot(t_check, y_actual, 'k--s', label="Actual y(t)", linewidth=1.2, markersize=3, zorder=10)
    plt.xlabel("t")
    plt.ylabel("State Value")
    plt.title("Numerical Methods Comparison for x(t) and y(t)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Phase plot: x(t) vs y(t)
    plt.figure(figsize=(10, 6))
    for m in active_methods:
        style, color = styles[m]
        plt.plot(series_x[m], series_y[m], style, label=method_labels[m], color=color, linewidth=0.8, markersize=3)
    if show_actual_in_plots:
        plt.plot(x_actual, y_actual, 'k--o', label="Actual", linewidth=0.8, markersize=4, zorder=10)
    plt.xlabel("x(t)")
    plt.ylabel("y(t)")
    plt.title("Phase Plot: x(t) vs y(t)")
    plt.legend()
    plt.grid(True)
    plt.show()

    if show_actual_in_plots:
        plt.figure(figsize=(10, 6))
        for m in active_methods:
            errs_x = [pct_err(x_actual[i], series_x[m][idx]) for i, idx in enumerate(indices)]
            marker = styles[m][0][0]
            color = styles[m][1]
            plt.scatter(t_check, errs_x, marker=marker, label=method_labels[m], color=color, zorder=5)
        plt.xlabel("t")
        plt.ylabel("Percent Relative Error (%)")
        plt.title("x(t) Percent Relative Errors")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

        plt.figure(figsize=(10, 6))
        for m in active_methods:
            errs_y = [pct_err(y_actual[i], series_y[m][idx]) for i, idx in enumerate(indices)]
            marker = styles[m][0][0]
            color = styles[m][1]
            plt.scatter(t_check, errs_y, marker=marker, label=method_labels[m], color=color, zorder=5)
        plt.xlabel("t")
        plt.ylabel("Percent Relative Error (%)")
        plt.title("y(t) Percent Relative Errors")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
else:
    print("matplotlib not installed; skipping plots.")
