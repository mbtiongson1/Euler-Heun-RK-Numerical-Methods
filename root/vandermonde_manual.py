"""Wrapper entrypoint for manual-point Vandermonde fitting.

Prefer: python -m numerical_methods.fitting.vandermonde_manual
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.fitting.vandermonde_manual")

