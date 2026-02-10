# utils.py

def print_table(headers, rows):
    """
    Print a table nicely formatted for console output.
    Floats are truncated to 6 decimal places.
    """
    print(" | ".join(headers))
    print("-" * (len(headers) * 12))
    for row in rows:
        formatted_row = []
        for val in row:
            if isinstance(val, float):
                formatted_row.append(f"{val:.6f}")
            else:
                formatted_row.append(str(val))
        print(" | ".join(formatted_row))


def print_table_csv(headers, rows):
    """
    Print a table in CSV format (comma-separated),
    ready for Excel copy-paste.
    Floats formatted to 6 decimals, percentages preserved as strings.
    """
    print(",".join(headers))
    for row in rows:
        formatted_row = []
        for val in row:
            if isinstance(val, str) and val.endswith('%'):
                formatted_row.append(val)
            elif isinstance(val, float):
                formatted_row.append(f"{val:.6f}")
            else:
                formatted_row.append(str(val))
        print(",".join(formatted_row))


def plot_polynomial(xs, ys, coeffs, label="Polynomial Fit", title="Polynomial Approximation"):
    """Plot the original data points and fitted polynomial curve.

    xs: array-like, x-coordinates of data points
    ys: array-like, y-coordinates of data points
    coeffs: array-like, polynomial coefficients (a_0, a_1, ..., a_p)
            where polynomial is a_0 + a_1*x + a_2*x^2 + ... + a_p*x^p
    label: label for the curve in the legend
    title: title of the plot
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        raise ImportError("matplotlib and numpy required for plotting.") from e

    xs_np = np.array(xs, dtype=float)
    ys_np = np.array(ys, dtype=float)
    coeffs_np = np.array(coeffs, dtype=float)

    # Create dense x grid for smooth polynomial curve
    x_min, x_max = xs_np.min(), xs_np.max()
    x_dense = np.linspace(x_min, x_max, 400)

    # Evaluate polynomial at dense points: sum(c_j * x^j)
    p_dense = np.polyval(coeffs_np[::-1], x_dense)  # polyval expects descending order

    plt.figure(figsize=(10, 6))
    plt.scatter(xs_np, ys_np, color='blue', s=50, label='Data Points', zorder=5)
    plt.plot(x_dense, p_dense, color='red', linewidth=2, label=label)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


def plot_polynomials_compare(xs, ys, coeffs_list, labels, title="Polynomial Comparison"):
    """Plot multiple fitted polynomials on the same axes for comparison.

    xs: array-like, x-coordinates of data points
    ys: array-like, y-coordinates of data points
    coeffs_list: list of coefficient arrays (each array is a_0..a_p for one fit)
    labels: list of labels for each polynomial curve
    title: title of the plot
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        raise ImportError("matplotlib and numpy required for plotting.") from e

    xs_np = np.array(xs, dtype=float)
    ys_np = np.array(ys, dtype=float)

    # Create dense x grid for smooth curves
    x_min, x_max = xs_np.min(), xs_np.max()
    x_dense = np.linspace(x_min, x_max, 400)

    plt.figure(figsize=(12, 7))
    plt.scatter(xs_np, ys_np, color='black', s=60, label='Data Points', zorder=10)

    colors = ['red', 'green', 'blue', 'orange', 'purple']
    for i, (coeffs, label) in enumerate(zip(coeffs_list, labels)):
        coeffs_np = np.array(coeffs, dtype=float)
        p_dense = np.polyval(coeffs_np[::-1], x_dense)
        color = colors[i % len(colors)]
        plt.plot(x_dense, p_dense, color=color, linewidth=2, label=label)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()