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

CMD_LINE = '[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] CollectGcBiasMetrics I=[BamDir]/[SeqID].analysisReady.bam O=[BamDir]/[SeqID].gc.bias.metrics.txt S=[BamDir]/[SeqID].gc.bias.summary.txt CHART=[BamDir]/[SeqID].gc.bias.metrics.pdf R=[ReferenceFasta] MINIMUM_GC=[min_gc] MAXIMUM_GC=[max_gc] WINDOW_SIZE=[window_size]'
ALL_KEYS = ['SeqID', 'BamDir', 'ReferenceFasta', 'java_bin', 'picard_jar', 'gc_threads', 'xmx_mb', 'min_gc', 'max_gc', 'window_size', 'gc_bias_metrics_txt', 'gc_bias_summary_txt', 'gc_bias_chart_pdf']
REQUIRED = {'SeqID': True, 'BamDir': True, 'ReferenceFasta': True, 'gc_bias_metrics_txt': False, 'gc_bias_summary_txt': False, 'gc_bias_chart_pdf': False}
DEFAULTS = {'java_bin': 'java', 'picard_jar': '/storage/apps/bin/picard-3.1.0.jar', 'gc_threads': '14', 'xmx_mb': '16384', 'min_gc': '0', 'max_gc': '100', 'window_size': '100', 'gc_bias_metrics_txt': '[BamDir]/[SeqID].gc.bias.metrics.txt', 'gc_bias_summary_txt': '[BamDir]/[SeqID].gc.bias.summary.txt', 'gc_bias_chart_pdf': '[BamDir]/[SeqID].gc.bias.metrics.pdf'}
OUTPUTS_META = {'gc_bias_metrics_txt': {'default': '[BamDir]/[SeqID].gc.bias.metrics.txt'}, 'gc_bias_summary_txt': {'default': '[BamDir]/[SeqID].gc.bias.summary.txt'}, 'gc_bias_chart_pdf': {'default': '[BamDir]/[SeqID].gc.bias.metrics.pdf'}}

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
        description="Runner: picard 3.1.0 (gc_bias_metrics)",
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
