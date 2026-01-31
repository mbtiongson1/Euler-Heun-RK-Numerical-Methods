import math

def f(x,y):
    return math.cos(x) - x*y

#initial conditions
x0 = 0
y0 = 1

#step size and range
h = 0.5
xn = 2

#actual solution to calculate error
y_actual = [1.0, 1.32326, 1.20079, 0.74740, 0.24241]