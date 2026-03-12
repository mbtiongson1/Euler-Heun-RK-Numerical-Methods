"""Wrapper entrypoint for running the fitting workflow driver.

Prefer: python -m numerical_methods.main.function
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.main.function")

