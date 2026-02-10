import math

def f(x,y):
    #return math.cos(x) - x*y for Probset 1
    return math.cos(x)-x*math.sin(x)

#initial conditions
x0 = 0
#y0 = 1 for Probset 1
y0 = 1

#step size and range
#h = 0.5
#xn = 2 for Probset 1
#xn = 3
m = 0 #change this for m-refinement
h=0.3/(2**m)

xn=1.5 #limit of x. By default, starts at 0.

#actual solution to calculate error

def y_actual_func(x):
    return x*math.cos(x) + 1

def generate_actual_solution(x0, xn, h, y_func):
    x = x0
    values = []
    while x <= xn + 1e-9:  # small tolerance for floating point
        values.append(y_func(x))
        x += h
    return values

y_actual = generate_actual_solution(x0, xn, h, y_actual_func)
#y_actual = [1.0, 1.32326, 1.20079, 0.74740, 0.24241] for Probset 1