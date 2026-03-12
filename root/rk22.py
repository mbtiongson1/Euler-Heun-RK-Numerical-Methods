"""Wrapper entrypoint for RK2.2 / Ralston's method (single IVP).

Prefer: python -m numerical_methods.methods.rk22
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.methods.rk22")

