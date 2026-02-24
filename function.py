"""Build Vandermonde, Lagrange, and Least-Squares polynomial approximations."""

from FD import compute_euler, compute_heun, compute_rk22, compute_rk4
from vandermonde import build_vandermonde, solve_vandermonde
from lagrange import solve_lagrange
from leastsquares import solve_least_squares
from ivp import method, p, xn, x0, y_actual, ls_methods
from utils import plot_polynomial, plot_polynomials_compare


METHOD_COMPUTE_MAP = {
    "euler": compute_euler,
    "heun": compute_heun,
    "rk22": compute_rk22,
    "rk4": compute_rk4,
}


def _poly_str(coeffs):
    terms = []
    for j, c in enumerate(coeffs):
        coeff_str = f"{c:.6f}"
        if j == 0:
            terms.append(coeff_str)
        else:
            terms.append(f"{coeff_str}*x^{j}")
    return " + ".join(terms)


def _eval_poly(coeffs, x):
    total = 0.0
    for j, c in enumerate(coeffs):
        total += c * (x ** j)
    return total


def _pct_err(actual, approx):
    if actual == 0:
        return 0.0
    return abs((actual - approx) / actual * 100)


def _print_ls_section(method_name, result):
    print("\n" + "=" * 60)
    print(f"LEAST-SQUARES ({method_name.upper()})")
    print("=" * 60)
    print(f"n (data points) = {result['n']}")
    print(f"p (polynomial degree) = {result['p']}")
    print("Candidate basis functions:")
    for j in range(result["p"] + 1):
        if j == 0:
            print("  phi_0(x) = 1")
        elif j == 1:
            print("  phi_1(x) = x")
        else:
            print(f"  phi_{j}(x) = x^{j}")
    print("Candidate model: y_hat(x) = a0 + a1*x + a2*x^2 + ... + ap*x^p")
    print("Loss: J(a) = sum_i (y_i - y_hat(x_i))^2 = || y - V a ||^2")

    print("\nNormal-equation matrix A = V^T V:")
    for row in result["A"]:
        print("  [" + ", ".join(f"{val:.6f}" for val in row) + "]")

    print("\nRight-hand side b = V^T y:")
    print("  [" + ", ".join(f"{val:.6f}" for val in result["b"]) + "]")

    print("\nLeast-squares Coefficients:")
    for idx, c in enumerate(result["coeffs"], start=1):
        print(f"a{idx} = {c:.6f}")

    print("\nLeast-squares Polynomial:")
    print(f"f(x) = {_poly_str(result['coeffs'])}")
    print(f"SSE = {result['sse']:.10f}")
    print(f"MSE = {result['mse']:.10f}")


def main():
    if method not in METHOD_COMPUTE_MAP:
        raise ValueError("Unknown method. Choose euler, heun, rk22, or rk4")

    xs, ys = METHOD_COMPUTE_MAP[method]()

    # Use only the first p+1 data points for interpolation-based methods
    xs_sampled, ys_sampled = xs[:p + 1], ys[:p + 1]

    has_actual = (
        isinstance(y_actual, list)
        and len(y_actual) > 1
        and all(isinstance(v, (int, float)) for v in y_actual)
    )
    if has_actual:
        # Compute x coordinates for y_actual reference points
        actual_spacing = (xn - x0) / (len(y_actual) - 1)
        x_actual = [x0 + i * actual_spacing for i in range(len(y_actual))]
    else:
        x_actual = None

    V = build_vandermonde(xs_sampled, p)

    print(f"Method: {method.upper()}, data points: {len(xs_sampled)}, polynomial degree p={p}")
    print("\nVandermonde Matrix (all rows, coefficients only):")
    for i in range(len(V)):
        x_val = V[i][0][0]
        y_val = ys_sampled[i]
        coeffs_only = [f"{val:.4f}" for _, val in V[i]]
        print(f"  [{', '.join(coeffs_only)}]  (x={x_val:.4f}, y={y_val:.4f})")
    print("Vandermonde shape:", (len(V), p + 1))

    coeffs = solve_vandermonde(V, ys_sampled)
    print("\nCoefficients:")
    for idx, c in enumerate(coeffs, start=1):
        print(f"a{idx} = {c:.6f}")
    print("\nVandermonde Polynomial:")
    print(f"f(x) = {_poly_str(coeffs)}")

    lagrange_coeffs = solve_lagrange(xs_sampled, ys_sampled, p=p)
    print("\n" + "=" * 60)
    print("LAGRANGE INTERPOLATION")
    print("=" * 60)
    print("\nIn terms of (x):")
    for idx, c in enumerate(lagrange_coeffs, start=1):
        print(f"b{idx} = {c:.6f}")
    print("\nLagrange Polynomial:")
    print(f"f(x) = {_poly_str(lagrange_coeffs)}")

    print("\n" + "=" * 60)
    print("COMPARISON PLOT")
    print("=" * 60)
    plot_polynomials_compare(
        xs_sampled,
        ys_sampled,
        [coeffs, lagrange_coeffs],
        [f"Vandermonde (degree {p})", f"Lagrange (degree {p})"],
        title="Vandermonde vs Lagrange Polynomial Approximation",
        x_end=xn,
        x_actual=x_actual if has_actual else None,
        y_actual=y_actual if has_actual else None,
    )
    plot_polynomial(
        xs_sampled,
        ys_sampled,
        coeffs,
        label="Vandermonde Fit",
        title=f"Vandermonde Polynomial Fit (p={p})",
        x_end=xn,
        x_actual=x_actual if has_actual else None,
        y_actual=y_actual if has_actual else None,
    )
    plot_polynomial(
        xs_sampled,
        ys_sampled,
        lagrange_coeffs,
        label="Lagrange Interpolation",
        title=f"Lagrange Interpolation (degree p={p})",
        x_end=xn,
        x_actual=x_actual if has_actual else None,
        y_actual=y_actual if has_actual else None,
    )

    # Least-squares comparison across multiple methods (configured in ivp.py)
    least_squares_results = []
    for ls_method in ls_methods:
        if ls_method not in METHOD_COMPUTE_MAP:
            continue
        xs_ls, ys_ls = METHOD_COMPUTE_MAP[ls_method]()
        try:
            ls_result = solve_least_squares(xs_ls, ys_ls, p)
            least_squares_results.append(
                {"method": ls_method, "xs": xs_ls, "ys": ys_ls, "result": ls_result}
            )
            _print_ls_section(ls_method, ls_result)
        except ValueError as e:
            print(f"\nSkipping least-squares for {ls_method.upper()}: {e}")

    if least_squares_results:
        coeffs_list = [entry["result"]["coeffs"] for entry in least_squares_results]
        labels = [
            f"Least-Squares {entry['method'].upper()} (degree {p})"
            for entry in least_squares_results
        ]
        xs_ref = least_squares_results[0]["xs"]
        ys_ref = least_squares_results[0]["ys"]
        plot_polynomials_compare(
            xs_ref,
            ys_ref,
            coeffs_list,
            labels,
            title="Least-Squares Fit Comparison Across Methods",
            x_end=xn,
            x_actual=x_actual if has_actual else None,
            y_actual=y_actual if has_actual else None,
            show_data_points=False,
        )

        # Error graph: least-squares fit vs y_actual at actual checkpoints
        if has_actual:
            try:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(10, 6))
                markers = ['o', 's', '^', 'd', 'p', '*']
                colors = ['blue', 'green', 'orange', 'red', 'purple', 'brown']

                for i, entry in enumerate(least_squares_results):
                    coeffs_ls = entry["result"]["coeffs"]
                    errors = []
                    for x_i, y_act in zip(x_actual, y_actual):
                        y_fit = _eval_poly(coeffs_ls, x_i)
                        errors.append(_pct_err(y_act, y_fit))

                    plt.scatter(
                        x_actual,
                        errors,
                        marker=markers[i % len(markers)],
                        label=f"LS {entry['method'].upper()}",
                        color=colors[i % len(colors)],
                        zorder=5,
                    )

                plt.xlabel("x")
                plt.ylabel("Percent Relative Error (%)")
                plt.title("Least-Squares Error vs Actual Values")
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.show()
            except Exception:
                print("matplotlib not installed; skipping least-squares error plot.")

    # Build least-squares error rows at y_actual checkpoints for CSV comparison
    ls_error_tables = []
    if has_actual:
        for entry in least_squares_results:
            method_name = entry["method"]
            coeffs_ls = entry["result"]["coeffs"]
            rows = []
            for i, (x_i, y_act) in enumerate(zip(x_actual, y_actual)):
                y_fit = _eval_poly(coeffs_ls, x_i)
                err_pct = _pct_err(y_act, y_fit)
                rows.append((i, x_i, y_act, y_fit, err_pct))
            ls_error_tables.append({"method": method_name, "rows": rows})

    _export_to_csv(
        xs_sampled,
        ys_sampled,
        V,
        coeffs,
        lagrange_coeffs,
        least_squares_results,
        ls_methods,
        ls_error_tables,
    )


def _export_to_csv(
    xs,
    ys,
    V,
    vand_coeffs,
    lag_coeffs,
    least_squares_results,
    ls_methods,
    ls_error_tables,
):
    """Export Vandermonde, Lagrange, and Least-Squares results to output.csv."""
    import csv

    V_values = [[val[1] for val in row] for row in V]

    with open("output.csv", "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow(["VANDERMONDE, LAGRANGE, AND LEAST-SQUARES APPROXIMATION"])
        writer.writerow([])
        writer.writerow(["Approximation methods used"])
        writer.writerow(["Interpolation method from ivp.py", method])
        writer.writerow(["Least-squares comparison methods", ", ".join(ls_methods)])
        writer.writerow(["Polynomial degree p", p])
        writer.writerow([])

        writer.writerow(["Data Points"])
        writer.writerow(["Index", "x", "y"])
        for i, (x, y) in enumerate(zip(xs, ys)):
            writer.writerow([i, f"{x:.6f}", f"{y:.6f}"])
        writer.writerow([])

        writer.writerow(["Vandermonde Matrix (V)"])
        writer.writerow(["Index", "x", "y"] + [f"x^{j}" for j in range(len(V_values[0]))])
        for i, (x, y) in enumerate(zip(xs, ys)):
            writer.writerow([i, f"{x:.6f}", f"{y:.6f}"] + [f"{val:.6f}" for val in V_values[i]])
        writer.writerow([])

        writer.writerow(["Vandermonde Coefficients"])
        writer.writerow(["Index", "Coefficient", "Value"])
        for idx, c in enumerate(vand_coeffs, start=1):
            writer.writerow([idx, f"a{idx}", f"{c:.6f}"])
        writer.writerow([])

        writer.writerow(["Lagrange Coefficients"])
        writer.writerow(["Index", "Coefficient", "Value"])
        for idx, c in enumerate(lag_coeffs, start=1):
            writer.writerow([idx, f"b{idx}", f"{c:.6f}"])
        writer.writerow([])

        writer.writerow(["Least-Squares Results"])
        if not least_squares_results:
            writer.writerow(["No least-squares results generated (likely p >= n)."])
        else:
            for entry in least_squares_results:
                method_name = entry["method"]
                result = entry["result"]
                writer.writerow([f"Method: {method_name.upper()}"])
                writer.writerow(["n", result["n"]])
                writer.writerow(["p", result["p"]])
                writer.writerow(["SSE", f"{result['sse']:.10f}"])
                writer.writerow(["MSE", f"{result['mse']:.10f}"])
                writer.writerow([])

                writer.writerow([f"A = V^T V ({method_name.upper()})"])
                for row in result["A"]:
                    writer.writerow([f"{val:.6f}" for val in row])
                writer.writerow([])

                writer.writerow([f"b = V^T y ({method_name.upper()})"])
                writer.writerow([f"{val:.6f}" for val in result["b"]])
                writer.writerow([])

                writer.writerow([f"Least-squares Coefficients ({method_name.upper()})"])
                writer.writerow(["Index", "Coefficient", "Value"])
                for idx, c in enumerate(result["coeffs"], start=1):
                    writer.writerow([idx, f"a{idx}", f"{c:.6f}"])
                writer.writerow([])

        writer.writerow(["Least-Squares Percent Relative Error vs y_actual"])
        if not ls_error_tables:
            writer.writerow(["No least-squares error table available."])
        else:
            # Combined comparison table: one row per x_actual, one error column per LS method
            methods_in_table = [entry["method"] for entry in ls_error_tables]
            header = ["Index", "x", "y_actual"]
            for m in methods_in_table:
                header.append(f"y_fit ({m.upper()})")
                header.append(f"Percent Relative Error ({m.upper()})")
            writer.writerow(header)

            n_rows = len(ls_error_tables[0]["rows"])
            for i in range(n_rows):
                base_row = ls_error_tables[0]["rows"][i]
                out_row = [base_row[0], f"{base_row[1]:.6f}", f"{base_row[2]:.6f}"]
                for method_table in ls_error_tables:
                    _, _, _, y_fit, err_pct = method_table["rows"][i]
                    out_row.append(f"{y_fit:.6f}")
                    out_row.append(f"{err_pct:.6f}%")
                writer.writerow(out_row)
            writer.writerow([])

    print("\n" + "=" * 60)
    print("EXPORT COMPLETE")
    print("=" * 60)
    print("Results exported to output.csv")


if __name__ == "__main__":
    main()
