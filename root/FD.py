"""Wrapper entrypoint for running finite differences.

Prefer: python -m numerical_methods.fd.FD
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.fd.FD")

