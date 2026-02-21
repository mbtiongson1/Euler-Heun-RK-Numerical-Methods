import math

def f(x,y):
    return -1*y**2-math.exp(-x)+x

#initial conditions
x0 = 0
y0 = 1
xn = 2.0 #limit of x. By default, starts at 0.

#Step Size, Number of steps, refinement parameters
m = 1 #change this for m-refinement
h=0.5/(2**m) #m-refinement
# h = 0.015625

#n = 9 #change this for n-refinement
#h=(xn/(n-1)) #n-refinement


# Numerical method and polynomial degree for approximation
method = 'rk4'  # Choose: euler, heun, rk22, rk4
p = 8  # Polynomial degree for Vandermonde and Lagrange fits

#actual solution to calculate error

def y_actual_func(x):
    return 0 #currently there is no actual function

def generate_actual_solution(x0, xn, h, y_func):
    x = x0
    values = []
    while x <= xn + 1e-9:  # small tolerance for floating point
        values.append(y_func(x))
        x += h
    return values

#y_actual = generate_actual_solution(x0, xn, h, y_actual_func)


def validate_data_points(num_points, polynomial_degree):
	num_intervals = num_points - 1  # Number of intervals between points
	
	if num_intervals % polynomial_degree == 0:
		# (num_points - 1) is a multiple of p
		step = num_intervals // polynomial_degree
		sampled_indices = list(range(0, num_points, step))[:polynomial_degree + 1]
		print(f"\nData Point Sampling:")
		print(f"  Total data points: {num_points} (indices 0 to {num_points-1})")
		print(f"  Polynomial degree: {polynomial_degree}")
		print(f"  Sampling every {step}-th point to get {len(sampled_indices)} points")
		print(f"  Sampled indices: {sampled_indices}")
		return sampled_indices, len(sampled_indices)
	else:
		# (num_points - 1) is not a multiple of p
		error_msg = (
			f"\nERROR: Invalid configuration!\n"
			f"  Number of intervals (n-1={num_intervals}) must be a multiple of p={polynomial_degree}\n"
			f"  Possible solutions:\n"
			f"    - Change p to a divisor of {num_intervals}\n"
			f"    - Change n to satisfy: (n-1) % p == 0\n"
			f"\n  Current options for p with n={num_points}:\n"
		)
		for p_candidate in range(num_intervals, -1, -1):
			if p_candidate == 0:
				continue
			if num_intervals % p_candidate == 0:
				error_msg += f"    - p={p_candidate} (step={num_intervals // p_candidate})\n"
		raise ValueError(error_msg)

