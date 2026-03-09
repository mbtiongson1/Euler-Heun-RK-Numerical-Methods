from ivpsystems2nd import f2, x0, y0, dy0, h, xn
from matrix import solve_linear_system
from utils import print_table


def linearize_f2(x):
    """
    Infer linear form coefficients for f2 at fixed x:
        f2(x, y, dy) = alpha(x)*y + beta(x)*dy + gamma(x)
    """
    gamma = f2(x, 0.0, 0.0)
    alpha = f2(x, 1.0, 0.0) - gamma
    beta = f2(x, 0.0, 1.0) - gamma

    # Consistency check to ensure f2 is linear in y and dy.
    test_points = [(2.0, -1.0), (-0.75, 1.25), (0.3, 0.4)]
    tol = 1e-9
    for y_test, dy_test in test_points:
        lhs = f2(x, y_test, dy_test)
        rhs = alpha * y_test + beta * dy_test + gamma
        if abs(lhs - rhs) > tol:
            raise ValueError(
                "f2 must be linear in (y, dy) for matrix-based finite differences.\n"
                f"Failed at x={x}, y={y_test}, dy={dy_test}: f2={lhs}, linear model={rhs}."
            )
    return alpha, beta, gamma


def build_fd_system():
    steps_exact = (xn - x0) / h
    n_steps = int(round(steps_exact))
    if abs(steps_exact - n_steps) > 1e-12:
        raise ValueError(
            "Inconsistent grid: h must divide (xn - x0) exactly. "
            f"Got (xn - x0)/h = {steps_exact}."
        )

    if n_steps == 0:
        return n_steps, [], []

    # Unknown vector u = [y1, y2, ..., yN]^T where N = n_steps.
    a = [[0.0] * n_steps for _ in range(n_steps)]
    b = [0.0] * n_steps

    # Discretized IC for derivative at x0 (forward difference):
    # (y1 - y0) / h = dy0  ->  y1 = y0 + h*dy0
    a[0][0] = 1.0
    b[0] = y0 + h * dy0

    # Discretized GE for interior nodes i = 1..N-1:
    # (y_{i+1} - 2y_i + y_{i-1}) / h^2 = alpha_i*y_i + beta_i*(y_{i+1}-y_{i-1})/(2h) + gamma_i
    # Rearranged:
    # (1 + beta_i*h/2)*y_{i-1} + (-2 - alpha_i*h^2)*y_i + (1 - beta_i*h/2)*y_{i+1} = gamma_i*h^2
    for i in range(1, n_steps):
        x_i = x0 + i * h
        alpha_i, beta_i, gamma_i = linearize_f2(x_i)
        a_left = 1.0 + 0.5 * beta_i * h
        a_center = -2.0 - alpha_i * h * h
        a_right = 1.0 - 0.5 * beta_i * h
        rhs = gamma_i * h * h

        row = i

        # y_{i-1}
        if i - 1 == 0:
            rhs -= a_left * y0
        else:
            a[row][i - 2] = a_left

        # y_i and y_{i+1}
        a[row][i - 1] = a_center
        a[row][i] = a_right
        b[row] = rhs

    return n_steps, a, b


def estimate_derivative(y_vals):
    n_steps = len(y_vals) - 1
    if n_steps == 0:
        return [dy0]

    dy_vals = []
    for i in range(n_steps + 1):
        if i == 0:
            dy_i = dy0
        elif i == n_steps:
            dy_i = (y_vals[i] - y_vals[i - 1]) / h
        else:
            dy_i = (y_vals[i + 1] - y_vals[i - 1]) / (2.0 * h)
        dy_vals.append(dy_i)
    return dy_vals


def print_matrix(a, b):
    print("\n3) Matrix System A*u = b (u = [y1, y2, ..., yN]^T)")
    print("A =")
    for row in a:
        print("  [" + ", ".join(f"{val:.6f}" for val in row) + "]")
    print("b = [" + ", ".join(f"{val:.6f}" for val in b) + "]")


def main():
    n_steps, a, b = build_fd_system()

    print("Finite Difference Method for Second-Order IVP")
    print(f"Problem setup: y'' = f2(x, y, y'), x in [{x0}, {xn}], h = {h}")
    print(f"Initial conditions: y({x0}) = {y0}, y'({x0}) = {dy0}")

    print("\n1) Discretized Governing Equation (interior node i = 1..N-1)")
    print("(y_{i+1} - 2y_i + y_{i-1})/h^2 = alpha_i*y_i + beta_i*(y_{i+1}-y_{i-1})/(2h) + gamma_i")
    print("(1 + beta_i*h/2)*y_{i-1} + (-2 - alpha_i*h^2)*y_i + (1 - beta_i*h/2)*y_{i+1} = gamma_i*h^2")

    print("\n2) Discretized Initial Conditions")
    print(f"y_0 = {y0}")
    print(f"(y_1 - y_0)/h = dy0  ->  y_1 = y_0 + h*dy0 = {y0 + h * dy0:.6f}")

    if n_steps == 0:
        print("\nNo unknowns to solve because xn == x0.")
        rows = [(0, x0, y0, dy0)]
        print("\n5) Estimated Table")
        print_table(["i", "x_i", "y_i", "y'_i"], rows)
        return

    print_matrix(a, b)

    u = solve_linear_system(a, b)
    y_vals = [y0] + u

    print("\n4) Unknown Values")
    for i, y_i in enumerate(y_vals[1:], start=1):
        print(f"y_{i} = {y_i:.10f}")

    x_vals = [x0 + i * h for i in range(n_steps + 1)]
    dy_vals = estimate_derivative(y_vals)

    rows = []
    for i, (x_i, y_i, dy_i) in enumerate(zip(x_vals, y_vals, dy_vals)):
        rows.append((i, x_i, y_i, dy_i))

    print("\n5) Estimated Table from Unknown Values")
    print_table(["i", "x_i", "y_i", "y'_i"], rows)


if __name__ == "__main__":
    main()
