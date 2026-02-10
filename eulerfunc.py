# euler.py

from ivp import (
    f, x0, y0, xn,
    generate_actual_solution,
    y_actual_func
)

def euler(h):
    x = x0
    y = y0

    xs = [x]
    ys = [y]

    # generate exact solution on this grid
    y_actual = generate_actual_solution(x0, xn, h, y_actual_func)
    es = [0.0]

    n = 0
    while x < xn and n < len(y_actual) - 1:
        y = y + h * f(x, y)
        x = x + h
        n += 1

        et = abs((y_actual[n] - y) / y_actual[n]) * 100

        xs.append(x)
        ys.append(y)
        es.append(et)

    return xs, ys, es