import math

def f(x,y):
    return math.cos(x)-x*math.sin(y)

#initial conditions
x0 = 0
y0 = 1

#step size and range

m = 1 #change this for m-refinement
#n = 5 #change this for n-refinement
h=0.3/(2**m) #m-refinement
#h=(xn/n) #n-refinement

xn=1.500000 #limit of x. By default, starts at 0.

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

#y_actual = generate_actual_solution(x0, xn, h, y_actual_func)