# Repo Architecture Context

This repository is organized as an importable Python package, with optional wrapper scripts under `root/`.

## High-Level Layout

- `numerical_methods/`: Canonical source code (import from here).
- `root/`: Optional wrapper scripts (for `python root/solver.py` instead of `python -m ...`).
- `docs/`: Static HTML index and documentation assets.
- `out/`: Generated outputs (CSVs, figures, reports). This folder is ignored by git.

## Package Responsibilities

- `numerical_methods/problems/`
  - Problem definitions and configuration values you edit (for example: `f`, initial conditions, step sizes, and any `*_actual` arrays).

- `numerical_methods/methods/`
  - Single-method scripts for solving IVPs/systems/shooting using a specific numerical method.

- `numerical_methods/main/`
  - Main “driver” scripts that compare multiple methods or orchestrate a workflow:
    - `solver.py`: single-IVP multi-method comparison
    - `systems.py`: system-IVP multi-method comparison
    - `shooting.py`: shooting-problem comparison
    - `function.py`: fitting workflow (Vandermonde/Lagrange/least-squares)

- `numerical_methods/fitting/`
  - Polynomial fitting building blocks and standalone fitting scripts.

- `numerical_methods/fd/`
  - Finite difference workflows and related scripts.

- `numerical_methods/experiments/`
  - Scratch/analysis scripts that are useful for exploration but are not core library code.

## How To Run

Preferred (module runs):
- `python -m numerical_methods.main.solver`
- `python -m numerical_methods.main.function`

Optional wrappers:
- `python root/solver.py`
- `python root/function.py`

## Outputs

Scripts that write files should write to `out/` (for example: `out/csv/`) via `numerical_methods/paths.py` so the repo root stays clean.
