"""1-D finite element problem definition.

Edit this file to define the boundary value problem solved by
``python -m numerical_methods.fea.fea1d``.
"""

problem_name = "Sample 1-D Poisson BVP"

# Domain [x0, xn]
x0 = 0.0
xn = 1.0


# Coefficient/source functions for: -(k(x)u'(x))' + c(x)u(x) = s(x)
def k(x):
    return 1.0


def c(x):
    return 0.0


def s(x):
    return 1.0


# Boundary conditions.
# Dirichlet format: {"type": "dirichlet", "value": ...}
# Neumann format:   {"type": "neumann", "value": ...}
# For Neumann, "value" is the boundary term based on k(x)u'(x).
left_bc = {"type": "dirichlet", "value": 0.0}
right_bc = {"type": "dirichlet", "value": 0.0}


# Mesh controls.
# mesh_mode choices:
#   "m"        -> h_used = base_h / (2 ** m)
#   "h"        -> use explicit h
#   "elements" -> use explicit num_elements
#   "manual"   -> use manual_nodes exactly
mesh_mode = "elements"
base_h = 0.25
m = 0
h = 0.25
num_elements = 4
manual_nodes = None


# Solver/output controls
quadrature_order = 2
print_level = "stage"  # "stage", "verbose", or "final"
export_csv = True
plot_result = True


def exact_solution(x):
    return 0.5 * x * (1.0 - x)
