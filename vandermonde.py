"""Vandermonde matrix builder and solver utilities."""


def build_vandermonde(xs, p):
    """Build Vandermonde matrix with tuples (x, x^j) for display."""
    return [[(x, x ** j) for j in range(p + 1)] for x in xs]


def solve_vandermonde(V, y):
    """Solve V a = y in least-squares sense and return coefficients."""
    try:
        import numpy as np
    except Exception as e:
        raise ImportError("NumPy required") from e

    # Extract polynomial values from tuples if present
    V_values = []
    for row in V:
        if row and isinstance(row[0], tuple):
            V_values.append([val[1] for val in row])
        else:
            V_values.append(row)

    V_np = np.array(V_values, dtype=float)
    y_np = np.array(y, dtype=float)
    coeffs, *_ = np.linalg.lstsq(V_np, y_np, rcond=None)
    return coeffs


if __name__ == '__main__':
	"""Standalone mode: compute Vandermonde polynomial fit from current IVP config."""
	from FD import compute_heun, compute_euler, compute_rk22, compute_rk4
	from ivp import method, p, validate_data_points

	# Compute x and y values using the configured method
	if method == 'euler':
		xs, ys = compute_euler()
	elif method == 'heun':
		xs, ys = compute_heun()
	elif method == 'rk22':
		xs, ys = compute_rk22()
	elif method == 'rk4':
		xs, ys = compute_rk4()
	else:
		raise ValueError(f'Unknown method: {method}')

	# Validate data points for polynomial degree
	sampled_indices, num_sampled = validate_data_points(len(xs), p)
	
	# Sample the data points
	xs_sampled = [xs[i] for i in sampled_indices]
	ys_sampled = [ys[i] for i in sampled_indices]

	# Build matrix and solve
	V = build_vandermonde(xs_sampled, p)
	coeffs = solve_vandermonde(V, ys_sampled)

	# Display results
	print("="*60)
	print("VANDERMONDE POLYNOMIAL FIT")
	print("="*60)
	print(f"Method: {method.upper()}, Degree: {p}")
	print(f"Data points: {len(xs_sampled)}")

	# Print coefficients as a1, a2, ...
	print("\nCoefficients:")
	for idx, c in enumerate(coeffs, start=1):
		print(f"a{idx} = {c:.6f}")

	# Build polynomial string
	terms = []
	for j, c in enumerate(coeffs):
		coeff_str = f"{c:.6f}"
		if j == 0:
			terms.append(f"{coeff_str}")
		else:
			terms.append(f"{coeff_str}*x^{j}")
	poly_str = " + ".join(terms)

	print(f"\nVandermonde Polynomial:")
	print(f"f(x) = {poly_str}")
