import csv

from ivpshooting import f, g, t0 as x0, x0 as y0, y0 as z0, h, tn as xn, x_actual as y_actual, y_actual as z_actual
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
num_steps = int(round((xn - x0) / h))
x_vals = [x0]

# ---------------------------
# Containers (all numerical method steps)
# ---------------------------
euler_y = [y0]
euler_z = [z0]

heun_y = [y0]
heun_z = [z0]

ral_y = [y0]
ral_z = [z0]

rk3_y = [y0]
rk3_z = [z0]

rk4_y = [y0]
rk4_z = [z0]

# ---------------------------
# Main loop
# ---------------------------
y_eu = y_he = y_ra = y_rk3 = y_rk4 = y0
z_eu = z_he = z_ra = z_rk3 = z_rk4 = z0

for step in range(num_steps):
    x = x0 + step * h
    x_next = x0 + (step + 1) * h

    # Euler
    kx = h * f(y_eu, z_eu, x)
    ky = h * g(y_eu, z_eu, x)
    y_eu = y_eu + kx
    z_eu = z_eu + ky

    # Heun
    yp = y_he + h * f(y_he, z_he, x)
    zp = z_he + h * g(y_he, z_he, x)
    y_he = y_he + (h / 2) * (f(y_he, z_he, x) + f(yp, zp, x_next))
    z_he = z_he + (h / 2) * (g(y_he, z_he, x) + g(yp, zp, x_next))

    # Ralston (RK2.2)
    k1x = h * f(y_ra, z_ra, x)
    k1y = h * g(y_ra, z_ra, x)
    k2x = h * f(y_ra + 3 * k1x / 4, z_ra + 3 * k1y / 4, x + 3 * h / 4)
    k2y = h * g(y_ra + 3 * k1x / 4, z_ra + 3 * k1y / 4, x + 3 * h / 4)
    y_ra = y_ra + (1 / 3) * k1x + (2 / 3) * k2x
    z_ra = z_ra + (1 / 3) * k1y + (2 / 3) * k2y

    # RK3 (Kutta's 3rd order)
    k1x = h * f(y_rk3, z_rk3, x)
    k1y = h * g(y_rk3, z_rk3, x)
    k2x = h * f(y_rk3 + k1x / 2, z_rk3 + k1y / 2, x + h / 2)
    k2y = h * g(y_rk3 + k1x / 2, z_rk3 + k1y / 2, x + h / 2)
    k3x = h * f(y_rk3 - k1x + 2 * k2x, z_rk3 - k1y + 2 * k2y, x + h)
    k3y = h * g(y_rk3 - k1x + 2 * k2x, z_rk3 - k1y + 2 * k2y, x + h)
    y_rk3 = y_rk3 + (1 / 6) * (k1x + 4 * k2x + k3x)
    z_rk3 = z_rk3 + (1 / 6) * (k1y + 4 * k2y + k3y)

    # RK4
    k1x = h * f(y_rk4, z_rk4, x)
    k1y = h * g(y_rk4, z_rk4, x)
    k2x = h * f(y_rk4 + k1x / 2, z_rk4 + k1y / 2, x + h / 2)
    k2y = h * g(y_rk4 + k1x / 2, z_rk4 + k1y / 2, x + h / 2)
    k3x = h * f(y_rk4 + k2x / 2, z_rk4 + k2y / 2, x + h / 2)
    k3y = h * g(y_rk4 + k2x / 2, z_rk4 + k2y / 2, x + h / 2)
    k4x = h * f(y_rk4 + k3x, z_rk4 + k3y, x_next)
    k4y = h * g(y_rk4 + k3x, z_rk4 + k3y, x_next)
    y_rk4 = y_rk4 + (1 / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    z_rk4 = z_rk4 + (1 / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)

    x_vals.append(x_next)
    euler_y.append(y_eu)
    euler_z.append(z_eu)
    heun_y.append(y_he)
    heun_z.append(z_he)
    ral_y.append(y_ra)
    ral_z.append(z_ra)
    rk3_y.append(y_rk3)
    rk3_z.append(z_rk3)
    rk4_y.append(y_rk4)
    rk4_z.append(z_rk4)


def pct_err(actual, approx):
    if actual == 0:
        return 0.0
    return abs((actual - approx) / actual * 100)


series_y = {
    "euler": euler_y,
    "heun": heun_y,
    "ral": ral_y,
    "rk3": rk3_y,
    "rk4": rk4_y,
}

series_z = {
    "euler": euler_z,
    "heun": heun_z,
    "ral": ral_z,
    "rk3": rk3_z,
    "rk4": rk4_z,
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
    isinstance(y_actual, list)
    and isinstance(z_actual, list)
    and len(y_actual) > 1
    and len(y_actual) == len(z_actual)
)

if has_actual:
    n_actual = len(y_actual)
    actual_spacing = (xn - x0) / (n_actual - 1)
    actual_stride = max(1, round(actual_spacing / h))
    x_check = [x0 + i * actual_spacing for i in range(n_actual)]
    indices = [i * actual_stride for i in range(n_actual)]
else:
    x_check = x_vals
    indices = list(range(len(x_vals)))

# Step: validate whether actual arrays are truly usable for plotting.
show_actual_in_plots = (
    has_actual
    and all(isinstance(v, (int, float)) for v in y_actual)
    and all(isinstance(v, (int, float)) for v in z_actual)
    and len(x_check) == len(y_actual) == len(z_actual)
)


# ---------------------------
# Tables
# ---------------------------
if has_actual:
    rows_y = []
    rows_z = []
    for i, idx in enumerate(indices):
        row_y = [i, x_check[i], y_actual[i]]
        row_z = [i, x_check[i], z_actual[i]]
        for m in active_methods:
            approx_y = series_y[m][idx]
            approx_z = series_z[m][idx]
            row_y.extend([approx_y, f"{pct_err(y_actual[i], approx_y):.5f}%"])
            row_z.extend([approx_z, f"{pct_err(z_actual[i], approx_z):.5f}%"])
        rows_y.append(tuple(row_y))
        rows_z.append(tuple(row_z))

    headers_y = ["n", "x", "y (actual)"]
    headers_z = ["n", "x", "z (actual)"]
    for m in active_methods:
        label = method_labels[m]
        headers_y.extend([f"{label} (y)", f"{label} (error)"])
        headers_z.extend([f"{label} (z)", f"{label} (error)"])

    print("Comparison of Numerical Methods for y(x)")
    print_table(headers_y, rows_y)
    print_table_csv(headers_y, rows_y)

    print("\nComparison of Numerical Methods for z(x)")
    print_table(headers_z, rows_z)
    print_table_csv(headers_z, rows_z)

    with open("output_systems_y.csv", mode="w", newline="") as fy:
        writer = csv.writer(fy)
        writer.writerow(headers_y)
        writer.writerows(rows_y)

    with open("output_systems_z.csv", mode="w", newline="") as fz:
        writer = csv.writer(fz)
        writer.writerow(headers_z)
        writer.writerows(rows_z)

    print("\nCSV files 'output_systems_y.csv' and 'output_systems_z.csv' successfully created!")
else:
    headers = ["n", "x"]
    for m in active_methods:
        label = method_labels[m]
        headers.extend([f"{label} (y)", f"{label} (z)"])

    rows = []
    for i, x in enumerate(x_vals):
        row = [i, x]
        for m in active_methods:
            row.extend([series_y[m][i], series_z[m][i]])
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

    y_palette = ['#1f77b4', '#2a9df4', '#1769aa', '#3b82f6', '#0ea5e9']
    z_palette = ['#d62728', '#ef4444', '#b91c1c', '#f97316', '#fb7185']
    label_stride = max(1, len(x_vals) // 12)

    plt.figure(figsize=(12, 7))
    for idx, m in enumerate(active_methods):
        y_color = y_palette[idx % len(y_palette)]
        z_color = z_palette[idx % len(z_palette)]

        plt.plot(
            x_vals, series_y[m], '-o',
            label=f"{method_labels[m]} y(x)",
            color=y_color, linewidth=1.0, markersize=3
        )
        plt.plot(
            x_vals, series_z[m], '--s',
            label=f"{method_labels[m]} z(x,y,z)",
            color=z_color, linewidth=1.0, markersize=3
        )

        for i in range(0, len(x_vals), label_stride):
            plt.annotate(
                f"{series_y[m][i]:.2f}",
                xy=(x_vals[i], series_y[m][i]),
                xytext=(0, 6),
                textcoords='offset points',
                fontsize=7,
                color=y_color,
                ha='center'
            )
            plt.annotate(
                f"{series_z[m][i]:.2f}",
                xy=(x_vals[i], series_z[m][i]),
                xytext=(0, -10),
                textcoords='offset points',
                fontsize=7,
                color=z_color,
                ha='center'
            )

    if show_actual_in_plots:
        plt.plot(x_check, y_actual, 'k-o', label="Actual y(x)", linewidth=1.2, markersize=3, zorder=10)
        plt.plot(x_check, z_actual, 'k--s', label="Actual z(x)", linewidth=1.2, markersize=3, zorder=10)
    plt.xlabel("x")
    plt.ylabel("State Value")
    plt.title("Numerical Methods Comparison for y(x) and z(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Phase plot: y(x) vs z(x)
    plt.figure(figsize=(10, 6))
    for m in active_methods:
        style, color = styles[m]
        plt.plot(series_y[m], series_z[m], style, label=method_labels[m], color=color, linewidth=0.8, markersize=3)
    if show_actual_in_plots:
        plt.plot(y_actual, z_actual, 'k--o', label="Actual", linewidth=0.8, markersize=4, zorder=10)
    plt.xlabel("y(x)")
    plt.ylabel("z(x)")
    plt.title("Phase Plot: y(x) vs z(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

    if show_actual_in_plots:
        plt.figure(figsize=(10, 6))
        for m in active_methods:
            errs_y = [pct_err(y_actual[i], series_y[m][idx]) for i, idx in enumerate(indices)]
            marker = styles[m][0][0]
            color = styles[m][1]
            plt.scatter(x_check, errs_y, marker=marker, label=method_labels[m], color=color, zorder=5)
        plt.xlabel("x")
        plt.ylabel("Percent Relative Error (%)")
        plt.title("y(x) Percent Relative Errors")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

        plt.figure(figsize=(10, 6))
        for m in active_methods:
            errs_z = [pct_err(z_actual[i], series_z[m][idx]) for i, idx in enumerate(indices)]
            marker = styles[m][0][0]
            color = styles[m][1]
            plt.scatter(x_check, errs_z, marker=marker, label=method_labels[m], color=color, zorder=5)
        plt.xlabel("x")
        plt.ylabel("Percent Relative Error (%)")
        plt.title("z(x) Percent Relative Errors")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
else:
    print("matplotlib not installed; skipping plots.")
