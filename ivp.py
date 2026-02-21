import math

def f(x,y):
    return math.sqrt(x) * math.sin(2*x) - 5*y

#initial conditions
x0 = 0
y0 = 0
xn = 2.4 #limit of x. By default, starts at 0.

#Step Size, Number of steps, refinement parameters
m = 1 #change this for m-refinement
h=0.3/(2**m) #m-refinement
# h = 0.3

#n = 9 #change this for n-refinement
#h=(xn/(n-1)) #n-refinement


# Numerical method and polynomial degree for approximation
method = 'rk22'  # Choose: euler, heun, rk22, rk3, rk4
p = 10  # Polynomial degree for Vandermonde and Lagrange fits

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
y_actual = [0, 0.025915, 0.094698, 0.156909, 0.166067, 0.101289, -0.023772, -0.165376, -0.266513]


# def validate_data_points(num_points, polynomial_degree):
# 	num_intervals = num_points - 1
# 	if num_intervals % polynomial_degree == 0:
# 		step = num_intervals // polynomial_degree
# 		sampled_indices = list(range(0, num_points, step))[:polynomial_degree + 1]
# 		return sampled_indices, len(sampled_indices)
# 	else:
# 		error_msg = (
# 			f"\nERROR: Invalid configuration!\n"
# 			f"  Number of intervals (n-1={num_intervals}) must be a multiple of p={polynomial_degree}\n"
# 		)
# 		raise ValueError(error_msg)

