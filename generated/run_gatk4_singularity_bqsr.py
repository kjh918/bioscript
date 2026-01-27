#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict

TOKEN_PAT = re.compile(r"\[([A-Za-z_][A-Za-z0-9_]*)\]")

CMD_LINE = 'singularity exec -B [bind] [sif] gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --reference [ReferenceFasta] --output [qcResDir]/[SeqID].recal.table.txt --known-sites [KnownSnp] --known-sites [KnownIndel1] --known-sites [KnownIndel2] && singularity exec -B [bind] [sif] gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m" --input [RecalInputBAM] --bqsr-recal-file [qcResDir]/[SeqID].recal.table.txt --output [BamDir]/[SeqID].recal.bam'
ALL_KEYS = ['SeqID', 'RecalInputBAM', 'ReferenceFasta', 'BamDir', 'qcResDir', 'KnownSnp', 'KnownIndel1', 'KnownIndel2', 'bind', 'sif', 'gc_threads', 'xmx_mb', 'recal_table', 'recal_bam', 'recal_bai']
REQUIRED = {'SeqID': True, 'RecalInputBAM': True, 'ReferenceFasta': True, 'BamDir': True, 'qcResDir': True, 'KnownSnp': True, 'KnownIndel1': True, 'KnownIndel2': True, 'recal_table': False, 'recal_bam': False, 'recal_bai': False}
DEFAULTS = {'bind': '/storage,/data', 'sif': '/storage/images/gatk-4.4.0.0.sif', 'gc_threads': '14', 'xmx_mb': '16384', 'recal_table': '[qcResDir]/[SeqID].recal.table.txt', 'recal_bam': '[BamDir]/[SeqID].recal.bam', 'recal_bai': '[BamDir]/[SeqID].recal.bam.bai'}
OUTPUTS_META = {'recal_table': {'default': '[qcResDir]/[SeqID].recal.table.txt'}, 'recal_bam': {'default': '[BamDir]/[SeqID].recal.bam'}, 'recal_bai': {'default': '[BamDir]/[SeqID].recal.bam.bai'}}

def render(s: str, ctx: Dict[str, str]) -> str:
    def repl(m: re.Match) -> str:
        key = m.group(1)
        return (ctx.get(key) or "").strip()
    out = TOKEN_PAT.sub(repl, s)
    out = re.sub(r"\s+", " ", out).strip()
    return out

def finalize_outputs(ctx: Dict[str, str]) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    for name, meta in OUTPUTS_META.items():
        meta = meta or {}
        cur = (ctx.get(name) or "").strip()
        if cur:
            outputs[name] = cur
            continue

        dv = meta.get("default")
        if isinstance(dv, str) and dv.strip():
            val = render(dv, ctx)
            outputs[name] = val
            ctx[name] = val
            continue

        if bool(meta.get("required") is True):
            raise SystemExit(f"Missing required output --{name} (no default/template inferred)")
    return outputs

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Runner: gatk4 4.4.0.0 (singularity_bqsr)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    for k in ALL_KEYS:
        if REQUIRED.get(k, False) and k not in OUTPUTS_META:
            ap.add_argument(f"--{k}", required=True, default="")
        else:
            ap.add_argument(f"--{k}", required=False, default=DEFAULTS.get(k, ""))

    ap.add_argument("--cwd", default=".")
    ap.add_argument("--print-cmd", action="store_true")
    ap.add_argument("--print-outputs", action="store_true")
    ap.add_argument("--emit-outputs", default="", help="write outputs json to file path")

    args = ap.parse_args()
    ctx = {k: ("" if v is None else str(v)) for k, v in vars(args).items()
            if k not in ("cwd","print_cmd","print_outputs","emit_outputs")}

    for k, dv in DEFAULTS.items():
        if not ctx.get(k, "").strip():
            ctx[k] = str(dv)

    outputs = finalize_outputs(ctx)

    cmd = render(CMD_LINE, ctx)

    if args.print_cmd:
        print(cmd)
        return

    proc = subprocess.run(shlex.split(cmd), cwd=str(Path(args.cwd)))
    rc = proc.returncode

    if args.print_outputs:
        print(json.dumps(outputs, ensure_ascii=False))

    if args.emit_outputs:
        Path(args.emit_outputs).write_text(json.dumps(outputs, ensure_ascii=False, indent=2), encoding="utf-8")

    raise SystemExit(rc)

if __name__ == "__main__":
    main()
