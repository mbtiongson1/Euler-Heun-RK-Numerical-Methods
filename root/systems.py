"""Wrapper entrypoint for running the system-IVP comparison driver.

Prefer: python -m numerical_methods.main.systems
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.main.systems")

