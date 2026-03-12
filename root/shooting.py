"""Wrapper entrypoint for running the shooting comparison driver.

Prefer: python -m numerical_methods.main.shooting
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.main.shooting")

