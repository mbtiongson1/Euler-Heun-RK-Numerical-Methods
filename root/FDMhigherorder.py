"""Wrapper entrypoint for running higher-order finite difference / FDM workflow.

Prefer: python -m numerical_methods.fd.FDMhigherorder
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.fd.FDMhigherorder")

