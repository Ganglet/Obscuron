"""Loads YAML config from the project's config/ directory."""

from __future__ import annotations

from pathlib import Path

import yaml

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"


def load_yaml(name: str) -> dict:
    path = CONFIG_DIR / name
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def models_config() -> dict:
    return load_yaml("models.yaml")


def snapshots_config() -> dict:
    return load_yaml("snapshots.yaml")
