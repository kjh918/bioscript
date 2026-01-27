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

CMD_LINE = 'singularity exec -B [bind] [sif] fastp --thread [Threads] --in1 [RawFastqDir]/[SeqID]_R1.fastq.gz --in2 [RawFastqDir]/[SeqID]_R2.fastq.gz --out1 [out_read1] --out2 [out_read2] --json [json] --html [html] --trim_poly_g --detect_adapter_for_pe --adapter_sequence [adapter_sequence] --adapter_sequence_r2 [adapter_sequence_r2] --length_required [length_required] --average_qual [average_qual] --qualified_quality_phred [qualified_quality_phred] --trim_front1 [trim_front1] --trim_front2 [trim_front2]'
ALL_KEYS = ['SeqID', 'RawFastqDir', 'TrimFastqDir', 'qcResDir', 'Threads', 'sif', 'bind', 'length_required', 'average_qual', 'qualified_quality_phred', 'trim_front1', 'trim_front2', 'adapter_sequence', 'adapter_sequence_r2', 'out_read1', 'out_read2', 'json', 'html']
REQUIRED = {'SeqID': True, 'RawFastqDir': True, 'TrimFastqDir': True, 'qcResDir': True, 'out_read1': False, 'out_read2': False, 'json': False, 'html': False}
DEFAULTS = {'Threads': '8', 'sif': '/storage/images/fastp-0.23.4.sif', 'bind': '/storage,/data', 'length_required': '100', 'average_qual': '10', 'qualified_quality_phred': '15', 'trim_front1': '14', 'trim_front2': '14', 'adapter_sequence': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'adapter_sequence_r2': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'out_read1': '[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz', 'out_read2': '[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz', 'json': '[qcResDir]/[SeqID].fastp.json', 'html': '[qcResDir]/[SeqID].fastp.html'}
OUTPUTS_META = {'out_read1': {'default': '[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz'}, 'out_read2': {'default': '[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz'}, 'json': {'default': '[qcResDir]/[SeqID].fastp.json'}, 'html': {'default': '[qcResDir]/[SeqID].fastp.html'}}

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
        description="Runner: fastp 0.23.4 (PicoplexGold_singularity_pe_trim)",
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
