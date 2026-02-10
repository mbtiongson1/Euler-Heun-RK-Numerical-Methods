def lagrange_basis(xs, i):
    n = len(xs)
    # Start with polynomial 1
    poly = [1.0]

    for j in range(n):
        if j == i:
            continue
        denom = xs[i] - xs[j]
        # Multiply poly by (x - xs[j]) / denom
        new_poly = [0.0] * (len(poly) + 1)
        for k in range(len(poly)):
            new_poly[k + 1] += poly[k] / denom
            new_poly[k] -= poly[k] * xs[j] / denom
        poly = new_poly

    return poly


def solve_lagrange(xs, ys, p=None):
    try:
        import numpy as np
    except Exception:
        raise ImportError("NumPy is required for solve_lagrange.")

    n = len(xs)
    if p is None:
        p = n - 1

    if p > n:
        # If requested degree exceeds number of points, use max interpolation
        max_degree = n - 1
    else:
        # Use requested degree p
        max_degree = p

    # Compute all Lagrange basis polynomials
    basis_polys = []
    for i in range(n):
        basis = lagrange_basis(xs, i)
        # Pad to max_degree
        basis.extend([0.0] * (max_degree + 1 - len(basis)))
        basis_polys.append(basis)

    # Combine: P(x) = sum_i ys[i] * L_i(x)
    combined = np.zeros(max_degree + 1)
    for i in range(n):
        combined += ys[i] * np.array(basis_polys[i])

    return combined.tolist()


if __name__ == '__main__':
	from FD import compute_euler, compute_heun, compute_rk22, compute_rk4
	from ivp import method, p, validate_data_points

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

	# Validate data points for polynomial degree
	sampled_indices, num_sampled = validate_data_points(len(xs), p)
	
	# Sample the data points
	xs_sampled = [xs[i] for i in sampled_indices]
	ys_sampled = [ys[i] for i in sampled_indices]

	print(f"Lagrange Interpolation using {method.upper()}, polynomial degree p={p}")
	print(f"Data points: {len(xs_sampled)}")

	# Compute Lagrange interpolation with degree p
	lagrange_coeffs = solve_lagrange(xs_sampled, ys_sampled, p=p)

	# Print coefficients as b_i values
	print("\nCoefficients:")
	for idx, c in enumerate(lagrange_coeffs, start=1):
		print(f"b{idx} = {c:.6f}")

	# Build Lagrange polynomial string
	lag_terms = []
	for j, c in enumerate(lagrange_coeffs):
		coeff_str = f"{c:.6f}"
		if j == 0:
			lag_terms.append(f"{coeff_str}")
		else:
			lag_terms.append(f"{coeff_str}*x^{j}")
	lag_poly_str = " + ".join(lag_terms)

	print("\nLagrange Polynomial:")
	print(f"f(x) = {lag_poly_str}")
