from __future__ import annotations

import runpy
import sys
from pathlib import Path


def _repo_root() -> Path:
    # This file lives at <repo>/root/_root_bootstrap.py
    return Path(__file__).resolve().parent.parent


def _ensure_repo_root_on_syspath() -> None:
    repo_root = _repo_root()
    repo_root_str = str(repo_root)
    if repo_root_str not in sys.path:
        sys.path.insert(0, repo_root_str)


def run(module: str) -> None:
    _ensure_repo_root_on_syspath()
    runpy.run_module(module, run_name="__main__")
