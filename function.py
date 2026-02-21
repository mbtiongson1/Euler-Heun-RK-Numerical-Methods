"""Build Vandermonde matrix from numerical solution methods.

Usage: python function.py [method] [p]
  method: euler | heun | rk22 | rk4   (default: heun)
  p: polynomial order (degree), integer (default: 5)

This script imports the step generators from the solver modules
(`FD.py` provides compatible `compute_*` functions) and builds the
Vandermonde matrix V (N x (p+1)) with V[i,j] = x_i**j and the
corresponding right-hand side vector y. The matrix and vector are
printed (first rows) and the full shapes are displayed.
"""

from FD import compute_euler, compute_heun, compute_rk22, compute_rk4
from vandermonde import build_vandermonde, solve_vandermonde
from lagrange import solve_lagrange
from ivp import method, p, xn, x0, y_actual
from utils import plot_polynomial, plot_polynomials_compare


def print_matrix_sample(V, y, nmax=8):
	n = min(len(V), nmax)
	print("First %d rows of Vandermonde matrix (rows) and y:" % n)
	for i in range(n):
		print(i, ":", [f"{val:.6g}" for val in V[i]], " | y=", f"{y[i]:.6g}")


def main():
	# Import method and p from ivp.py (already imported at module level)

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

	# # Validate that there are enough data points for polynomial degree
	# sampled_indices, num_sampled = validate_data_points(len(xs), p)
	# xs_sampled = [xs[i] for i in sampled_indices]
	# ys_sampled = [ys[i] for i in sampled_indices]
	# Use only the first p+1 data points (enough for degree-p polynomial)
	xs_sampled, ys_sampled = xs[:p + 1], ys[:p + 1]

	# Compute x coordinates for y_actual reference points
	actual_spacing = xn / (len(y_actual) - 1)
	x_actual = [x0 + i * actual_spacing for i in range(len(y_actual))]

	V = build_vandermonde(xs_sampled, p)

	print(f"Method: {method.upper()}, data points: {len(xs_sampled)}, polynomial degree p={p}")
	print("\nVandermonde Matrix (all rows, coefficients only):")
	for i in range(len(V)):
		# V[i] is list of tuples (x, x^j)
		x_val = V[i][0][0]  # Get x from first tuple
		y_val = ys_sampled[i]  # Get corresponding y value
		coeffs_only = [f"{val:.4f}" for _, val in V[i]]
		print(f"  [{', '.join(coeffs_only)}]  (x={x_val:.4f}, y={y_val:.4f})")
	print("Vandermonde shape:", (len(V), p + 1))

	# Solve for polynomial coefficients using the shared vandermonde solver
	coeffs = solve_vandermonde(V, ys_sampled)

	# Print coefficients in user requested format: a1..an where a1 is coeff of x^0
	print("\nCoefficients:")
	for idx, c in enumerate(coeffs, start=1):
		print(f"a{idx} = {c:.6f}")

	# Build polynomial string f(x) = a1 + a2*x + a3*x^2 + ...
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


	# Compute Lagrange interpolation with degree p
	lagrange_coeffs = solve_lagrange(xs_sampled, ys_sampled, p=p)

	print("\n" + "="*60)
	print("LAGRANGE INTERPOLATION")
	print("="*60)

	# Print Lagrange coefficients as b_i values
	print("\nIn terms of (x):")
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
	

	# Compare both methods on same plot
	print("\n" + "="*60)
	print("COMPARISON PLOT")
	print("="*60)
	plot_polynomials_compare(xs_sampled, ys_sampled, [coeffs, lagrange_coeffs],
	                          [f"Vandermonde (degree {p})", f"Lagrange (degree {p})"],
	                          title="Vandermonde vs Lagrange Polynomial Approximation", x_end=xn,
	                          x_actual=x_actual, y_actual=y_actual)
	# Plot Vandermonde fit
	plot_polynomial(xs_sampled, ys_sampled, coeffs, label="Vandermonde Fit", title=f"Vandermonde Polynomial Fit (p={p})", x_end=xn,
	                x_actual=x_actual, y_actual=y_actual)
	# Plot Lagrange fit
	plot_polynomial(xs_sampled, ys_sampled, lagrange_coeffs, label="Lagrange Interpolation", title=f"Lagrange Interpolation (degree (p={p})", x_end=xn,
	                x_actual=x_actual, y_actual=y_actual)

	# Export results to CSV
	_export_to_csv(xs_sampled, ys_sampled, V, coeffs, lagrange_coeffs)


def _export_to_csv(xs, ys, V, vand_coeffs, lag_coeffs):
	"""Export data and results to output.csv for easier analysis."""
	import csv

	# Extract just polynomial values from tuples
	V_values = [[val[1] for val in row] for row in V]

	with open('output.csv', 'w', newline='') as f:
		writer = csv.writer(f)

		# Header
		writer.writerow(['VANDERMONDE AND LAGRANGE POLYNOMIAL APPROXIMATION'])
		writer.writerow([])

		# Data points section
		writer.writerow(['Data Points'])
		writer.writerow(['Index', 'x', 'y'])
		for i, (x, y) in enumerate(zip(xs, ys)):
			writer.writerow([i, f'{x:.6f}', f'{y:.6f}'])
		writer.writerow([])

		# Vandermonde matrix section
		writer.writerow(['Vandermonde Matrix (V)'])
		writer.writerow(['Index', 'x', 'y'] + [f'x^{j}' for j in range(len(V_values[0]))])
		for i, (x, y) in enumerate(zip(xs, ys)):
			writer.writerow([i, f'{x:.6f}', f'{y:.6f}'] + [f'{val:.6f}' for val in V_values[i]])
		writer.writerow([])

		# Vandermonde coefficients section
		writer.writerow(['Vandermonde Coefficients'])
		writer.writerow(['Index', 'Coefficient', 'Value'])
		for idx, c in enumerate(vand_coeffs, start=1):
			writer.writerow([idx, f'a{idx}', f'{c:.6f}'])
		writer.writerow([])

		# Lagrange coefficients section
		writer.writerow(['Lagrange Coefficients'])
		writer.writerow(['Index', 'Coefficient', 'Value'])
		for idx, c in enumerate(lag_coeffs, start=1):
			writer.writerow([idx, f'b{idx}', f'{c:.6f}'])

	print("\n" + "="*60)
	print("EXPORT COMPLETE")
	print("="*60)
	print("Results exported to output.csv")


if __name__ == '__main__':
	main()
