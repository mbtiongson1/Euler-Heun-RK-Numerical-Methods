#RK4 method

from ivp import f, x0, y0, h, xn#, y_actual
from utils import print_table

#RK4 method

from ivp import f, x0, y0, h, xn#, y_actual
from utils import print_table

# Calculate number of steps
num_steps = int(round((xn - x0) / h))

rows = []
rows.append((0, x0, 0, 0, 0, 0, y0, 0))

y = y0
for n in range(1, num_steps + 1):
	x = x0 + n * h
	# Use previous x for k1, k2, k3
	x_prev = x0 + (n - 1) * h
	
	k1 = h*f(x_prev, y)
	k2 = h*f(x_prev + (1/2)*h, y+ (1/2)*k1)
	k3 = h*f(x_prev + (1/2)*h, y+ (1/2)*k2)
	k4 = h*f(x, y + k3)
	y = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
	rows.append((n, x, k1, k2, k3, k4, y))

print("Runge-Kutta 4th Order Method (RK4)")
print_table(["n", "x", "k1", "k2", "k3", "k4", "y"], rows)

#print("Runge-Kutta 4th Order Method (RK4)")
#print_table(["n", "x", "k1", "k2", "k3", "k4", "y"], rows)