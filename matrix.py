from ivpsecond import x0, y0, dy0, h
from utils import print_table
import math


def solve_linear_system(a, b, pivot_tol=1e-14):
    """
    Solve A x = b using Gaussian elimination with partial pivoting.

    Args:
        a: Coefficient matrix as a list of lists (n x n).
        b: Right-hand side vector (length n).
        pivot_tol: Minimum allowed absolute pivot magnitude.

    Returns:
        x: Solution vector (length n).
    """
    n = len(a)
    if n == 0:
        return []
    if len(b) != n:
        raise ValueError("Dimension mismatch: len(b) must equal number of rows in a.")
    for row in a:
        if len(row) != n:
            raise ValueError("Matrix a must be square (n x n).")

    aug = [row[:] + [rhs] for row, rhs in zip(a, b)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < pivot_tol:
            raise ValueError("Singular matrix encountered while solving A x = b.")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]

        for row in range(col + 1, n):
            factor = aug[row][col] / aug[col][col]
            if factor == 0.0:
                continue
            for k in range(col, n + 1):
                aug[row][k] -= factor * aug[col][k]

    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        rhs = aug[row][n] - sum(aug[row][k] * x[k] for k in range(row + 1, n))
        x[row] = rhs / aug[row][row]
    return x


def matrix(a, b, pivot_tol=1e-14):
    """Backward-compatible alias to the linear solver."""
    return solve_linear_system(a, b, pivot_tol=pivot_tol)


def solution_rows_from_coefficients(coeffs):
    """
    Build (i, x_i, y_i, y'_i) rows from solved coefficients.
    coeffs is assumed to be [y_1, y_2, ..., y_N].
    """
    y_vals = [y0] + list(coeffs)
    x_vals = [x0 + i * h for i in range(len(y_vals))]

    dy_vals = []
    n_steps = len(y_vals) - 1
    for i in range(n_steps + 1):
        if i == 0:
            dy_i = dy0
        elif i == n_steps:
            dy_i = (y_vals[i] - y_vals[i - 1]) / h
        else:
            dy_i = (y_vals[i + 1] - y_vals[i - 1]) / (2.0 * h)
        dy_vals.append(dy_i)

    rows = []
    for i, (x_i, y_i, dy_i) in enumerate(zip(x_vals, y_vals, dy_vals)):
        rows.append((i, x_i, y_i, dy_i))
    return rows


def print_solution_table(coeffs):
    rows = solution_rows_from_coefficients(coeffs)
    print_table(["i", "x_i", "y_i", "y'_i"], rows)


def main():
    h = 0.2
    a1 = 50
    a2 = 25
    a3 = 10*math.exp(-h)

    a = [
        [4, -1, 0.0, 0.0, 0.0],
        [h-a1, a2, 0.0, 0.0, 0.0],
        [a2, 2*h-a1, a2, 0.0, 0.0],
        [0.0, a2, 3*h-a1, a2, 0.0],
        [0.0, 0.0, a2, 4*h-a1, a2],
    ]

    b = [4.8, 10*math.exp(-h)-50, 10*math.exp(-2*h), 10*math.exp(-3*h), 10*math.exp(-4*h)]
    x = solve_linear_system(a, b)
    print("Coefficients:", x)
    print_solution_table(x)


if __name__ == "__main__":
    main()
