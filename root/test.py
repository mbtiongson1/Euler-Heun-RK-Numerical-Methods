"""Wrapper entrypoint for experiments/test script.

Prefer: python -m numerical_methods.experiments.test
"""

from _root_bootstrap import run


if __name__ == "__main__":
    run("numerical_methods.experiments.test")

