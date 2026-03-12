# Plan Archive: Restructure Repo Into `numerical_methods/` Package + Root Wrappers

Date: 2026-03-13 (Asia/Manila)

Snapshot: 0e58f004d5c820f14399880a3147e58463e4099b

Branch: restructure-folders

---

# Restructure Repo Into `numerical_methods/` Package + Root Wrappers

## Summary
Create a clean, importable Python package `numerical_methods/` with purpose-based subfolders, keep the current “many scripts” UX via root wrapper scripts, standardize generated outputs into `out/`, and document the architecture in `context.md`. Work happens on a new branch `restructure-folders`, with the current plan archived in `plan.md` and indexed to the pre-change snapshot commit.

## Snapshot + Branching (Before Any Code Moves)
1. Verify working tree is clean (`git status -sb`).
2. Create and switch to a new branch: `git checkout -b restructure-folders`.
3. Capture snapshot commit hash (the plan’s reference point): `git rev-parse HEAD`.
4. Add `plan.md` at repo root containing:
- Title + date
- `Snapshot: <commit-hash-from-step-3>`
- `Branch: restructure-folders`
- A full copy of this plan content (so it always refers to the exact snapshot it was written against)

## Target Structure (Canonical Code Location)
```
.
├── numerical_methods/
│   ├── __init__.py
│   ├── paths.py
│   ├── utils.py
│   ├── matrix.py
│   ├── problems/
│   ├── methods/
│   ├── main/
│   ├── fitting/
│   ├── fd/
│   └── experiments/
├── docs/
│   └── index.html
├── out/
│   ├── csv/
│   ├── figures/
│   └── reports/
├── context.md
└── (root wrapper scripts)
```

## File Moves / Renames (Decision-Complete Mapping)
Move these files into the package (use `git mv` for all tracked moves/renames):

### 1. Problems (configs you edit)
- `ivp.py` -> `numerical_methods/problems/ivp.py`
- `ivpsystems.py` -> `numerical_methods/problems/ivpsystems.py`
- `ivpshooting.py` -> `numerical_methods/problems/ivpshooting.py`
- `ivphigherorder.py` -> `numerical_methods/problems/ivphigherorder.py`

### 2. Methods (single IVP / systems / shooting step scripts)
- `euler.py` -> `numerical_methods/methods/euler.py`
- `heun.py` -> `numerical_methods/methods/heun.py`
- `rk22.py` -> `numerical_methods/methods/rk22.py`
- `rk3.py` -> `numerical_methods/methods/rk3.py`
- `rk4.py` -> `numerical_methods/methods/rk4.py`
- `heunsystems.py` -> `numerical_methods/methods/heunsystems.py`
- `rk4systems.py` -> `numerical_methods/methods/rk4systems.py`
- `heunshooting.py` -> `numerical_methods/methods/heunshooting.py`
- `rk4shooting.py` -> `numerical_methods/methods/rk4shooting.py`

### 3. Main entrypoints (renamed from `main*` naming)
- `mainsolver.py` -> `numerical_methods/main/solver.py`
- `mainsystems.py` -> `numerical_methods/main/systems.py`
- `mainshooting.py` -> `numerical_methods/main/shooting.py`
- `function.py` -> `numerical_methods/main/function.py` (per request)

### 4. Fitting modules (library-ish + runnable modules)
- `vandermonde.py` -> `numerical_methods/fitting/vandermonde.py`
- `vandermonde_manual.py` -> `numerical_methods/fitting/vandermonde_manual.py`
- `lagrange.py` -> `numerical_methods/fitting/lagrange.py`
- `leastsquares.py` -> `numerical_methods/fitting/leastsquares.py`

### 5. Finite difference / higher-order scripts
- `FD.py` -> `numerical_methods/fd/FD.py`
- `FDMhigherorder.py` -> `numerical_methods/fd/FDMhigherorder.py`

### 6. Core shared modules at package root
- `utils.py` -> `numerical_methods/utils.py`
- `matrix.py` -> `numerical_methods/matrix.py` (package root, per request)

### 7. Experiments / scratch analyses
- `test.py` -> `numerical_methods/experiments/test.py`
- `eulerfunc.py` -> `numerical_methods/experiments/eulerfunc.py`

### 8. Docs
- `index.html` -> `docs/index.html` and update relative links so they still resolve from `docs/`:
- README link becomes `../README.md`
- image links become `../Figure 1 - ...png` and `../Figure 2 - ...png`

Add `__init__.py` files in:
- `numerical_methods/`
- `numerical_methods/problems/`
- `numerical_methods/methods/`
- `numerical_methods/main/`
- `numerical_methods/fitting/`
- `numerical_methods/fd/`
- `numerical_methods/experiments/`

## Import Policy (Package-Absolute Imports Only)
Update all internal imports to use `numerical_methods.*` (no more root-level imports):
- `from ivp import ...` -> `from numerical_methods.problems.ivp import ...`
- `from ivpsystems import ...` -> `from numerical_methods.problems.ivpsystems import ...`
- `from ivpshooting import ...` -> `from numerical_methods.problems.ivpshooting import ...`
- `from ivphigherorder import ...` -> `from numerical_methods.problems.ivphigherorder import ...`
- `from utils import ...` -> `from numerical_methods.utils import ...`
- `from matrix import ...` -> `from numerical_methods.matrix import ...`
- `from FD import compute_*` -> `from numerical_methods.fd.FD import compute_*`
- `from vandermonde import ...` -> `from numerical_methods.fitting.vandermonde import ...`
- `from lagrange import ...` -> `from numerical_methods.fitting.lagrange import ...`
- `from leastsquares import ...` -> `from numerical_methods.fitting.leastsquares import ...`

## Output Standardization (`out/`)
1. Add `numerical_methods/paths.py` (uses `pathlib.Path`) that:
- Resolves repo root from the file location (not CWD)
- Exposes `OUT_DIR`, `CSV_DIR`, `FIGURES_DIR`, `REPORTS_DIR`
- Provides `csv_path(filename)` that ensures `out/csv/` exists before returning the path

2. Update all CSV writers to use `csv_path(...)`:
- `numerical_methods/main/solver.py`: `out/csv/output.csv`
- `numerical_methods/main/systems.py`: `out/csv/output_systems_x.csv`, `out/csv/output_systems_y.csv`
- `numerical_methods/main/shooting.py`: `out/csv/output_systems_y.csv`, `out/csv/output_systems_z.csv`
- `numerical_methods/fd/FD.py`: `out/csv/FD.csv`
- `numerical_methods/main/function.py`: change export target to `out/csv/output_fit.csv` (avoid clobbering `output.csv`)

3. Update `.gitignore`:
- Add `.venv/`
- Add `out/`

## Root Wrapper Scripts (Preserve Existing UX + New Names)
Create small root-level wrapper scripts that call into package modules via `runpy.run_module(..., run_name="__main__")`.

1. Canonical root entrypoints (new names)
- `solver.py` -> runs `numerical_methods.main.solver`
- `systems.py` -> runs `numerical_methods.main.systems`
- `shooting.py` -> runs `numerical_methods.main.shooting`
- `function.py` -> runs `numerical_methods.main.function`

2. Back-compat root entrypoints (old names retained)
- `mainsolver.py` -> runs `numerical_methods.main.solver`
- `mainsystems.py` -> runs `numerical_methods.main.systems`
- `mainshooting.py` -> runs `numerical_methods.main.shooting`

3. Keep root wrappers for the remaining runnable scripts (same filenames as before), delegating to their new module paths:
- `euler.py`, `heun.py`, `rk22.py`, `rk3.py`, `rk4.py`
- `heunsystems.py`, `rk4systems.py`
- `heunshooting.py`, `rk4shooting.py`
- `FD.py`, `FDMhigherorder.py`
- `vandermonde.py`, `vandermonde_manual.py`, `lagrange.py`, `leastsquares.py`
- `test.py`, `eulerfunc.py`

Do not keep root shims for library/config modules:
- No root `utils.py`, `matrix.py`, or `ivp*.py` after the move (they live under `numerical_methods/`).

## Documentation
1. Add `context.md` at repo root:
- Explain folder responsibilities (`problems/`, `methods/`, `main/`, `fitting/`, `fd/`, `experiments/`)
- State that `numerical_methods/` is the source of truth and root scripts are wrappers
- Show the preferred run style (`python -m ...`) vs wrapper runs (`python solver.py`, etc.)
- Mention outputs go to `out/`

2. Update `README.md`:
- Replace the old “Project Structure” tree with the new one
- Add a “How to Run” section showing:
- Preferred module runs: `python -m numerical_methods.main.solver` (etc.)
- Wrapper runs: `python solver.py`, `python mainsolver.py` (compat), `python euler.py` (etc.)
- Add a “Configuration” note explaining problem definitions moved to `numerical_methods/problems/*`

## Test Plan (Post-Change Verification)
Run from repo root:
- `python -c "import numerical_methods; import numerical_methods.matrix; import numerical_methods.utils"`
- `python -m compileall numerical_methods`
- `python -m numerical_methods.main.solver`
- `python -m numerical_methods.main.function`
- `python solver.py`
- `python mainsolver.py`
- `python FD.py rk4`
- Confirm CSVs are created under `out/csv/` only

## Assumptions / Defaults
- Dependencies remain as-is (no attempt to remove `numpy`, `matplotlib`, `pandas`).
- “Many scripts” is preserved via wrappers; preferred modern usage is `python -m ...`.
- `function.py` export filename becomes `output_fit.csv` to avoid collision with `solver.py`’s `output.csv`.

