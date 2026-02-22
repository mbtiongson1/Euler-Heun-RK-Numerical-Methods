"""Manual-point Vandermonde polynomial fit and plot against y_actual."""

from ivp import x0, xn, y_actual
from utils import plot_polynomial
from vandermonde import build_vandermonde, solve_vandermonde


if __name__ == '__main__':
    # Manual points for Vandermonde matrix construction (degree 10 polynomial)
    manual_points = [
        (0.0, 0.0),
        (0.3, 0.029901),
        (0.6, 0.096455),
        (0.75, 0.129568),
        (0.9, 0.1543211),
        (1.2, 0.160079),
        (1.5, 0.094608),
        (1.8, -0.027978),
        (2.1, -0.164812),
        (2.25, -0.22121),
        (2.4, -0.260755),
    ]

    xs_manual = [pt[0] for pt in manual_points]
    ys_manual = [pt[1] for pt in manual_points]
    p_manual = len(xs_manual) - 1

    # Build Vandermonde matrix and solve for coefficients
    V = build_vandermonde(xs_manual, p_manual)
    coeffs = solve_vandermonde(V, ys_manual)

    # Display results (same style as vandermonde.py)
    print("=" * 60)
    print("VANDERMONDE MANUAL-POINT POLYNOMIAL FIT")
    print("=" * 60)
    print(f"Degree: {p_manual}")
    print(f"Data points: {len(xs_manual)}")

    print("\nCoefficients:")
    for idx, c in enumerate(coeffs, start=1):
        print(f"a{idx} = {c:.6f}")

    terms = []
    for j, c in enumerate(coeffs):
        coeff_str = f"{c:.6f}"
        if j == 0:
            terms.append(f"{coeff_str}")
        else:
            terms.append(f"{coeff_str}*x^{j}")
    poly_str = " + ".join(terms)

    print("\nVandermonde Polynomial:")
    print(f"f(x) = {poly_str}")

    # Compute x-coordinates corresponding to y_actual checkpoints for plotting comparison
    actual_spacing = (xn - x0) / (len(y_actual) - 1)
    x_actual = [x0 + i * actual_spacing for i in range(len(y_actual))]

    # Plot manual-point polynomial fit and compare to actual values
    plot_polynomial(
        xs_manual,
        ys_manual,
        coeffs,
        label="Vandermonde Fit (Manual Points)",
        title=f"Manual Vandermonde Polynomial Fit (p={p_manual})",
        x_end=xn,
        x_actual=x_actual,
        y_actual=y_actual,
    )
