#!/usr/bin/env python
"""CLI: append one entry to docs/experiment_log.md.

Usage:
    uv run python scripts/log_experiment.py \\
        --title "full-panel snapshot differencing" \\
        --config "R207 panel (502 genomes), Pfam-35 GA / Pfam-37 net-new-family proxy" \\
        --result "1,341,100 proteins; 4,138 positive-proxy (0.309%)" \\
        --next "Report to Track 1; formal go/no-go still pending his sign-off"
"""

from __future__ import annotations

import argparse

from darkmatter.experiment_log import log_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--title", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--next", required=True, dest="next_step")
    args = parser.parse_args()

    log_experiment(args.title, args.config, args.result, args.next_step)
    print("appended entry to docs/experiment_log.md")


if __name__ == "__main__":
    main()
