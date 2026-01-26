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

CMD_LINE = 'singularity exec -B [bind] [sif] fastqc [fastqc_args] --threads [threads] --outdir [qcResDir] [RawFastqDir]/[SeqID]_R1.fastq.gz [RawFastqDir]/[SeqID]_R2.fastq.gz'
ALL_KEYS = ['SeqID', 'RawFastqDir', 'qcResDir', 'bind', 'sif', 'threads', 'fastqc_args', 'r1_html', 'r2_html', 'r1_zip', 'r2_zip', 'r1_dir', 'r2_dir']
REQUIRED = {'SeqID': True, 'RawFastqDir': True, 'qcResDir': True, 'r1_html': False, 'r2_html': False, 'r1_zip': False, 'r2_zip': False, 'r1_dir': False, 'r2_dir': False}
DEFAULTS = {'bind': '/storage,/data', 'sif': '/storage/images/fastqc-0.12.1.sif', 'threads': '8', 'fastqc_args': '--extract', 'r1_html': '[qcResDir]/[SeqID]_R1_fastqc.html', 'r2_html': '[qcResDir]/[SeqID]_R2_fastqc.html', 'r1_zip': '[qcResDir]/[SeqID]_R1_fastqc.zip', 'r2_zip': '[qcResDir]/[SeqID]_R2_fastqc.zip', 'r1_dir': '[qcResDir]/[SeqID]_R1_fastqc', 'r2_dir': '[qcResDir]/[SeqID]_R2_fastqc'}
OUTPUTS_META = {'r1_html': {'default': '[qcResDir]/[SeqID]_R1_fastqc.html'}, 'r2_html': {'default': '[qcResDir]/[SeqID]_R2_fastqc.html'}, 'r1_zip': {'default': '[qcResDir]/[SeqID]_R1_fastqc.zip'}, 'r2_zip': {'default': '[qcResDir]/[SeqID]_R2_fastqc.zip'}, 'r1_dir': {'default': '[qcResDir]/[SeqID]_R1_fastqc'}, 'r2_dir': {'default': '[qcResDir]/[SeqID]_R2_fastqc'}}

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
        description="Runner: fastqc 0.12.1 (singularity_pe_qc_extract)",
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
