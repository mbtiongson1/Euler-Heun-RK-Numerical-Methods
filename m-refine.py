# m-refine.py

from ivp import h, y_actual_func
from euler2 import euler

def m_refine(m):
    Xs = []
    Ys = []
    Es = []

    # run Euler on refined grids
    for k in range(m + 1):
        hk = h / (2 ** k)
        x, y, e = euler(hk)
        Xs.append(x)
        Ys.append(y)
        Es.append(e)

    # coarsest grid
    x_ref = Xs[0]

    # header
    header = f"{'x':<10}{'y_actual':<14}"
    for k in range(m + 1):
        header += f"y{k}".ljust(14)
    header += "e_m"
    print(header)
    print("-" * len(header))

    # rows
    for i, xval in enumerate(x_ref):
        row = f"{xval:<10.5f}{y_actual_func(xval):<14.5f}"

        for k in range(m + 1):
            stride = 2 ** k
            row += f"{Ys[k][i * stride]:<14.5f}"

        row += f"{Es[-1][i * (2 ** m)]:.5f}"
        print(row)


if __name__ == "__main__":
    m = 5
    m_refine(m)