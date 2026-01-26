from __future__ import annotations

import argparse
from pathlib import Path

from bioscript.generator import make_runners


def main() -> None:
    ap = argparse.ArgumentParser(prog="bio-script")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_make = sub.add_parser("make", help="Generate run_*.py and run_*.sh from YAML config")
    p_make.add_argument("-c", "--config", required=True, help="tool config YAML")
    p_make.add_argument("-o", "--outdir", default="", help="override output_dir in YAML (optional)")

    args = ap.parse_args()

    if args.cmd == "make":
        make_runners(
            Path(args.config).resolve(),
            outdir_override=(Path(args.outdir).resolve() if args.outdir else None),
        )
