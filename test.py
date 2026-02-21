"""Compare RK4 results for different refinement levels (n=9 and n=5)."""

import math
from ivp import f, x0, y0, xn

def compute_rk4(n_points):
	"""Compute RK4 solution with n refinement level."""
	h = xn / (n_points - 1)
	num_steps = n_points - 1
	
	y = y0
	results = [(x0, y0)]
	
	for n in range(1, num_steps + 1):
		x = x0 + n * h
		x_prev = x0 + (n - 1) * h
		
		k1 = h * f(x_prev, y)
		k2 = h * f(x_prev + (1/2)*h, y + (1/2)*k1)
		k3 = h * f(x_prev + (1/2)*h, y + (1/2)*k2)
		k4 = h * f(x, y + k3)
		y = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
		results.append((x, y))
	
	return results

def align_results(results1, results2, label1, label2):
	"""Align two sets of (x, y) results and display in table format.
	
	Args:
		results1: List of (x, y) tuples from first run
		results2: List of (x, y) tuples from second run
		label1: Label for first run
		label2: Label for second run
	"""
	# Create dictionaries for quick lookup
	dict1 = {round(x, 10): y for x, y in results1}
	dict2 = {round(x, 10): y for x, y in results2}
	
	# Collect all unique x values
	all_x = sorted(set(dict1.keys()) | set(dict2.keys()))
	
	# Print header
	print(f"\n{'n':>4} | {'x':>8} | {label1:>15} | {'n':>4} | {label2:>15}")
	print("-" * 60)
	
	# Print rows
	for x_val in all_x:
		idx1 = None
		idx2 = None
		y1_val = None
		y2_val = None
		
		# Find indices and values
		for i, (x, y) in enumerate(results1):
			if abs(x - x_val) < 1e-9:
				idx1 = i
				y1_val = y
				break
		
		for i, (x, y) in enumerate(results2):
			if abs(x - x_val) < 1e-9:
				idx2 = i
				y2_val = y
				break
		
		# Format output
		n1_str = str(idx1) if idx1 is not None else "-"
		n2_str = str(idx2) if idx2 is not None else "-"
		y1_str = f"{y1_val:.6f}" if y1_val is not None else "-"
		y2_str = f"{y2_val:.6f}" if y2_val is not None else "-"
		
		print(f"{n1_str:>4} | {x_val:>8.4f} | {y1_str:>15} | {n2_str:>4} | {y2_str:>15}")

# Compute RK4 for both refinement levels
print("="*60)
print("RK4 COMPARISON: m=0 vs m=1")
print("="*60)

results_n5 = compute_rk4(5)
results_n9 = compute_rk4(9)

print(f"\nRK4 with m=0 (5 points, h={xn/(5-1):.6f}):")
print(f"Points: {len(results_n5)}")
print(f"\nRK4 with m=1 (9 points, h={xn/(9-1):.6f}):")
print(f"Points: {len(results_n9)}")

# Align and display
align_results(results_n5, results_n9, "y (m=0)", "y (m=1)")

# Show differences where both have values
print("\n" + "="*60)
print("DIFFERENCES (where both have values)")
print("="*60)
dict9 = {round(x, 10): y for x, y in results_n9}
dict5 = {round(x, 10): y for x, y in results_n5}
common_x = sorted(set(dict9.keys()) & set(dict5.keys()))

print(f"{'x':>8} | {'y(m=0)':>15} | {'y(m=1)':>15} | {'|Δy|':>12} | {'Rel. Error %':>12}")
print("-" * 70)
for x_val in common_x:
	y9 = dict9[x_val]
	y5 = dict5[x_val]
	diff = abs(y5 - y9)
	rel_error = abs(diff / y5 * 100) if y5 != 0 else 0
	print(f"{x_val:>8.4f} | {y5:>15.6f} | {y9:>15.6f} | {diff:>12.6f} | {rel_error:>12.6f}")