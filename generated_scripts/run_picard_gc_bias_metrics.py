#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = '[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] CollectGcBiasMetrics I=[BamDir]/[SeqID].analysisReady.bam O=[BamDir]/[SeqID].gc.bias.metrics.txt S=[BamDir]/[SeqID].gc.bias.summary.txt CHART=[BamDir]/[SeqID].gc.bias.metrics.pdf R=[ReferenceFasta] MINIMUM_GC=[min_gc] MAXIMUM_GC=[max_gc] WINDOW_SIZE=[window_size]'
REQUIRED_KEYS = ['SeqID', 'BamDir', 'ReferenceFasta']
DEFAULTS = {'java_bin': 'java', 'picard_jar': '/storage/apps/bin/picard-3.1.0.jar', 'gc_threads': '14', 'xmx_mb': '16384', 'min_gc': '0', 'max_gc': '100', 'window_size': '100', 'gc_bias_metrics_txt': '[BamDir]/[SeqID].gc.bias.metrics.txt', 'gc_bias_summary_txt': '[BamDir]/[SeqID].gc.bias.summary.txt', 'gc_bias_chart_pdf': '[BamDir]/[SeqID].gc.bias.metrics.pdf'}
OUTPUT_KEYS = ['gc_bias_metrics_txt', 'gc_bias_summary_txt', 'gc_bias_chart_pdf']

def render(s, ctx):
    def repl(m):
        key = m.group(1)
        return str(ctx.get(key, m.group(0))).strip()
    res = s
    for _ in range(3):
        if "[" not in res: break
        res = TOKEN_PAT.sub(repl, res)
    return " ".join(res.split())

def main():
    parser = argparse.ArgumentParser()
    for k in ['BamDir', 'ReferenceFasta', 'SeqID', 'gc_bias_chart_pdf', 'gc_bias_metrics_txt', 'gc_bias_summary_txt', 'gc_threads', 'java_bin', 'max_gc', 'min_gc', 'picard_jar', 'window_size', 'xmx_mb']:
        parser.add_argument(f"--{k}", default=DEFAULTS.get(k, ""))
    parser.add_argument("--cwd", default=".")
    parser.add_argument("--emit-outputs", default="")

    args = parser.parse_args()
    ctx = vars(args)

    for rk in REQUIRED_KEYS:
        if not ctx.get(rk):
            print(f"Error: Missing --{rk}")
            exit(1)

    # Output 경로 렌더링
    for k in OUTPUT_KEYS:
        ctx[k] = render(ctx[k], ctx)

    final_cmd = render(CMD_LINE, ctx)
    print(f"\n[RUNNING CMD]\n{final_cmd}\n")

    os.makedirs(args.cwd, exist_ok=True)
    proc = subprocess.run(final_cmd, shell=True, cwd=args.cwd)

    if args.emit_outputs:
        out_data = {k: ctx[k] for k in OUTPUT_KEYS}
        Path(args.emit_outputs).write_text(json.dumps(out_data, indent=2))

    exit(proc.returncode)

if __name__ == "__main__":
    main()
