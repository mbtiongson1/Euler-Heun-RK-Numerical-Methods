"""Wrapper entrypoint for running the single-IVP comparison driver.

Prefer: python -m numerical_methods.main.solver
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.main.solver")

