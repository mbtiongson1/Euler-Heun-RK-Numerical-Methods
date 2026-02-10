from ivp import f, x0, y0, h, xn #, y_actual
from utils import print_table

x = x0
y = y0

rows =[]
rows.append((0, x, y, 0))

n = 0
num_steps = int(round((xn - x0) / h))
for n in range(1, num_steps + 1):
    y = y + h * f(x, y)
    x = x0 + n * h  # Recalculate x from x0 to avoid accumulation errors
    # True percent relative error
    #et = abs((y_actual[n] - y) / y_actual[n] * 100)
    #et_str = f"{et:.5f}"+"%" #just to make it 5 decimal places
    #rows.append((n, x, y, y_actual[n], et_str))
    rows.append((n, x, y))

print("Euler's Method")
#print_table(["n", "x", "y", "y_actual", "e_t (%)"], rows)
print_table(["n", "x", "y"], rows)
answer = y
error = abs((y - 0.955900)/y)*100
print(f"Error is", f"{error:.5f} %")