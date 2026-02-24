# Euler, Heun, Runge-Kutta Numerical Methods with Polynomial Fitting

Numerical analysis toolkit for:
- Single first-order IVPs (`y' = f(x, y)`)
- Systems of first-order IVPs (`x' = f(x, y, t)`, `y' = g(x, y, t)`)
- Polynomial approximation and interpolation
- Method comparison, tabular output, CSV export, and plotting

## Current Capabilities

### Single-IVP Solvers
- Euler (`euler.py`)
- Heun / Predictor-Corrector (`heun.py`)
- Ralston RK2.2 (`rk22.py`)
- RK3 (`rk3.py`)
- RK4 (`rk4.py`)

### System-IVP Solvers
- Heun for systems (`heunsystems.py`)
- RK4 for systems (`rk4systems.py`)
- Multi-method system comparison runner (`mainsystems.py`)
- System configuration and optional exact/manual actual values (`ivpsystems.py`)

### Approximation / Analysis
- Vandermonde fitting (`vandermonde.py`)
- Manual Vandermonde workflow (`vandermonde_manual.py`)
- Least-squares comparison workflow (`leastsquares.py`)
- Lagrange interpolation (`lagrange.py`)
- Finite differences (`FD.py`)

### Orchestration / Utilities
- Main solver (`mainsolver.py`)
- Combined workflow (`function.py`)
- Shared formatting and plotting utilities (`utils.py`)
- Comparison checks (`test.py`)
- Core single-IVP configuration (`ivp.py`)

## Project Structure

```
FD.py
README.md
euler.py
eulerfunc.py
function.py
heun.py
heunsystems.py
index.html
ivp.py
ivpsystems.py
lagrange.py
leastsquares.py
mainsolver.py
mainsystems.py
rk22.py
rk3.py
rk4.py
rk4systems.py
test.py
utils.py
vandermonde.py
vandermonde_manual.py
```

## Quick Start

1. Configure your problem in `ivp.py` (single IVP) or `ivpsystems.py` (system IVP).
2. Run one of the drivers:
- `python function.py`
- `python mainsolver.py`
- `python mainsystems.py`
3. For individual method runs, execute the corresponding module directly (e.g., `python rk4.py`, `python rk4systems.py`).

## Updates Since `main` Baseline Commit

Baseline on `main`: `25868de` ("Resolve merge conflict: merge unit1-function-approximations into main")

All commits introduced on this branch lineage after that baseline:

1. `fa95001` Improve polynomial plot formatting and update numerical method configs
- Updated: `ivp.py`, `function.py`, `mainsolver.py`, `utils.py`
- Result: improved plotting format and method configuration handling

2. `45beb74` Added RK3 and actual-solution plot comparisons
- Added: `rk3.py`
- Updated: `function.py`, `ivp.py`, `lagrange.py`, `mainsolver.py`, `utils.py`, `vandermonde.py`
- Result: RK3 support and overlays/comparisons against actual solution points

3. `2661d6a` Added manual Vandermonde and least-squares workflows
- Added: `vandermonde_manual.py`, `leastsquares.py`
- Updated: `FD.py`, `function.py`, `ivp.py`, `mainsolver.py`, `utils.py`
- Result: expanded approximation workflows and supporting integration updates

4. `aa84789` Added system-IVP module set
- Added: `ivpsystems.py`, `heunsystems.py`, `rk4systems.py`, `mainsystems.py`
- Result: dedicated configuration and solvers for 2-equation IVP systems with comparison tooling

5. `5f87c16` Removed redundant matplotlib flag
- Updated: `mainsystems.py`
- Result: plotting path simplified by removing unnecessary matplotlib availability guard

6. `0199baa` Minor follow-up updates for system-IVP flow
- Updated: `ivp.py`, `ivpsystems.py`, `mainsystems.py`
- Result: manual override support for actual solution generation and small cleanup refinements

Net file-level diff vs `main` includes:
- Modified: `FD.py`, `function.py`, `ivp.py`, `lagrange.py`, `mainsolver.py`, `rk4.py`, `test.py`, `utils.py`, `vandermonde.py`
- Added: `heunsystems.py`, `ivpsystems.py`, `leastsquares.py`, `mainsystems.py`, `rk3.py`, `rk4systems.py`, `vandermonde_manual.py`

## Dependencies

- Python 3.7+
- NumPy
- Matplotlib

Install:

```bash
pip install numpy matplotlib
```
