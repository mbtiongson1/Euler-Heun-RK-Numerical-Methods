"""Least-squares polynomial approximation for cases where p < n."""

from FD import compute_euler, compute_heun, compute_rk22, compute_rk4
from ivp import method, p, x0, xn, y_actual
from utils import plot_polynomial


def _print_matrix(name, M):
    print(f"\n{name} =")
    for row in M:
        print("  [" + ", ".join(f"{val:.6f}" for val in row) + "]")


def solve_least_squares(xs, ys, p):
    """Solve polynomial least-squares fit and return matrices, coefficients, and fit metrics."""
    try:
        import numpy as np
    except Exception as e:
        raise ImportError("NumPy is required for least-squares fitting.") from e

    n = len(xs)
    if p >= n:
        raise ValueError(
            f"Least-squares setup requires p < n. Current values: p={p}, n={n}."
        )

    V = np.array([[x ** j for j in range(p + 1)] for x in xs], dtype=float)
    y = np.array(ys, dtype=float)
    A = V.T @ V
    b = V.T @ y
    coeffs = np.linalg.solve(A, b)
    y_hat = V @ coeffs
    residuals = y - y_hat
    sse = float(np.sum(residuals ** 2))
    mse = sse / n

    return {
        "V": V,
        "A": A,
        "b": b,
        "coeffs": coeffs,
        "y_hat": y_hat,
        "sse": sse,
        "mse": mse,
        "n": n,
        "p": p,
    }


def main():
    if method == 'euler':
        xs, ys = compute_euler()
    elif method == 'heun':
        xs, ys = compute_heun()
    elif method == 'rk22':
        xs, ys = compute_rk22()
    elif method == 'rk4':
        xs, ys = compute_rk4()
    else:
        raise ValueError('Unknown method. Choose euler, heun, rk22, or rk4')

    result = solve_least_squares(xs, ys, p)
    V = result["V"]
    A = result["A"]
    b = result["b"]
    coeffs = result["coeffs"]
    sse = result["sse"]
    mse = result["mse"]
    n = result["n"]

    print("=" * 60)
    print("LEAST-SQUARES POLYNOMIAL APPROXIMATION")
    print("=" * 60)
    print(f"Method: {method.upper()}")
    print(f"n (data points) = {n}")
    print(f"p (polynomial degree) = {p}")

    print("\nCandidate basis functions:")
    for j in range(p + 1):
        if j == 0:
            print("  phi_0(x) = 1")
        elif j == 1:
            print("  phi_1(x) = x")
        else:
            print(f"  phi_{j}(x) = x^{j}")

    print("\nCandidate model:")
    print("  y_hat(x) = a0 + a1*x + a2*x^2 + ... + ap*x^p")

    print("\nLoss functions:")
    print("  J(a) = sum_{i=1..n} (y_i - y_hat(x_i))^2")
    print("  J(a) = || y - V a ||^2")

    _print_matrix("V (design matrix)", V)
    _print_matrix("A = V^T V", A)
    print("\nb = V^T y =")
    print("  [" + ", ".join(f"{val:.6f}" for val in b) + "]")

    print("\nCoefficients:")
    for idx, c in enumerate(coeffs):
        print(f"a{idx + 1} = {c:.6f}")

    poly_terms = []
    for j, c in enumerate(coeffs):
        coeff_str = f"{c:.6f}"
        if j == 0:
            poly_terms.append(coeff_str)
        else:
            poly_terms.append(f"{coeff_str}*x^{j}")
    print("\nFinal polynomial:")
    print("f(x) = " + " + ".join(poly_terms))

    print("\nFit error summary:")
    print(f"SSE = {sse:.10f}")
    print(f"MSE = {mse:.10f}")

    # Compare fit against reference y_actual points from ivp.py
    actual_spacing = (xn - x0) / (len(y_actual) - 1)
    x_actual = [x0 + i * actual_spacing for i in range(len(y_actual))]
    plot_polynomial(
        xs,
        ys,
        coeffs,
        label=f"Least-Squares Fit (degree {p})",
        title=f"Least-Squares Polynomial Fit ({method.upper()}, p={p})",
        x_end=xn,
        x_actual=x_actual,
        y_actual=y_actual,
    )


if __name__ == '__main__':
    main()
