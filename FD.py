# Finite Difference Methods using y from the method used (euler, heun, rk22, rk4)
# to estimate numerical differentiation of y at each step.
# Forward Difference at first step, Backward Difference at last step,
# Central Difference for interior points.

from ivp import f, x0, y0, h, xn
from utils import print_table
import sys


def compute_euler():
	x = x0
	y = y0
	xs = [x]
	ys = [y]
	num_steps = int(round((xn - x0) / h))
	for n in range(1, num_steps + 1):
		y = y + h * f(x, y)
		x = x0 + n * h
		xs.append(x)
		ys.append(y)
	return xs, ys


def compute_heun():
	x = x0
	y = y0
	xs = [x]
	ys = [y]
	num_steps = int(round((xn - x0) / h))
	for n in range(1, num_steps + 1):
		yp = y + h * f(x, y)
		x_next = x0 + n * h
		yc_next = y + (h / 2) * (f(x, y) + f(x_next, yp))
		x = x_next
		y = yc_next
		xs.append(x)
		ys.append(y)
	return xs, ys


def compute_rk22():
	x = x0
	y = y0
	a1 = 1/3
	a2 = 2/3
	xs = [x]
	ys = [y]
	num_steps = int(round((xn - x0) / h))
	for n in range(1, num_steps + 1):
		k1 = h * f(x, y)
		k2 = h * f(x + (3/4) * h, y + (3/4) * k1)
		y = y + a1 * k1 + a2 * k2
		x = x0 + n * h
		xs.append(x)
		ys.append(y)
	return xs, ys


def compute_rk4():
	x = x0
	y = y0
	xs = [x]
	ys = [y]
	num_steps = int(round((xn - x0) / h))
	for n in range(1, num_steps + 1):
		k1 = h * f(x, y)
		k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
		k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
		k4 = h * f(x + h, y + k3)
		y = y + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
		x = x0 + n * h
		xs.append(x)
		ys.append(y)
	return xs, ys


def finite_differences(xs, ys):
	n_pts = len(ys)
	rows = []
	for i in range(n_pts):
		x = xs[i]
		y = ys[i]
		if i == 0:
			dydx = (ys[1] - ys[0]) / h
			method = 'forward'
		elif i == n_pts - 1:
			dydx = (ys[-1] - ys[-2]) / h
			method = 'backward'
		else:
			dydx = (ys[i + 1] - ys[i - 1]) / (2 * h)
			method = 'central'
		# second derivative (d2y/dx2) using central difference for interior points, forward/backward for endpoints
		if n_pts >= 3:
			if i == 0:
				d2ydx2 = (ys[0] - 2 * ys[1] + ys[2]) / (h * h)
			elif i == n_pts - 1:
				d2ydx2 = (ys[-3] - 2 * ys[-2] + ys[-1]) / (h * h)
			else:
				d2ydx2 = (ys[i + 1] - 2 * ys[i] + ys[i - 1]) / (h * h)
		else:
			# Not enough points for second derivative; set to None
			d2ydx2 = None

		rows.append((i, x, y, dydx, d2ydx2, method))
	return rows


def main():
	method = 'rk4' #edit this!!!
	if len(sys.argv) > 1:
		method = sys.argv[1].lower()

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

	rows = finite_differences(xs, ys)

	print(f"Finite Differences using {method.upper()} method")
	print_table(["n", "x", "y", "dy/dx", "d2y/dx2", "method"], rows)


if __name__ == '__main__':
	main()
