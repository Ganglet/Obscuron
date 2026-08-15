"""Append entries to docs/experiment_log.md — one per meaningful run, per
the format docs/reproducibility.md's "Experiment tracking" section
promises: commit hash, config, result summary, and what it implies next.
"""

from __future__ import annotations

import subprocess
from datetime import date
from pathlib import Path

LOG_PATH = Path(__file__).resolve().parents[2] / "docs" / "experiment_log.md"


def _git_commit() -> str:
    return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()


def log_experiment(title: str, config: str, result: str, next_step: str, log_path: Path = LOG_PATH) -> None:
    entry = (
        f"\n## {date.today().isoformat()} — {title}\n"
        f"- Commit: {_git_commit()}\n"
        f"- Config: {config}\n"
        f"- Result: {result}\n"
        f"- Next: {next_step}\n"
    )
    with log_path.open("a", encoding="utf-8") as f:
        f.write(entry)
