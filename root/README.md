# Wrapper Scripts (`root/`)

This folder contains optional wrapper scripts so you can run commands like `python root/solver.py` instead of `python -m ...`.

Why this exists:
- Keeps the repository root clean (no duplicated entrypoints next to the package).
- Preserves a simple `python <script>.py` workflow for coursework runs.

## Main Drivers

- `python root/solver.py`
- `python root/systems.py`
- `python root/shooting.py`
- `python root/function.py`

## Methods / Tools

- `python root/rk4.py` (and other method wrappers)
- `python root/FD.py rk4`

## Notes

- These wrappers modify `sys.path` to ensure the repo root is importable, so they work even if you run them from inside `root/`.
- Preferred style (no wrappers) remains:
  - `python -m numerical_methods.main.solver`
  - `python -m numerical_methods.methods.rk4`

