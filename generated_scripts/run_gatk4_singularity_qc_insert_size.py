#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'export LC_ALL=[lc_all] && singularity exec -B [bind] [sif] [java_bin] -XX:ParallelGCThreads=[Threads]    -Xmx[xmx_mb]m -jar [gatk_jar] CollectInsertSizeMetrics --INPUT [BamDir]/[SeqID].analysisReady.bam --OUTPUT [qcResDir]/[SeqID].insert.size.metrics.txt --Histogram_FILE [qcResDir]/[SeqID].insert.size.histogram.pdf --REFERENCE_SEQUENCE [ReferenceFasta]'
REQUIRED_KEYS = ['SeqID', 'BamDir', 'qcResDir', 'ReferenceFasta']
DEFAULTS = {'bind': '/storage,/data', 'sif': '/storage/images/gatk-4.4.0.0.sif', 'Threads': '14', 'xmx_mb': '16384', 'gatk_jar': '/gatk/gatk-package-4.4.0.0-local.jar', 'java_bin': 'java', 'lc_all': 'en_US.UTF-8', 'insert_size_metrics_txt': '[qcResDir]/[SeqID].insert.size.metrics.txt', 'insert_size_hist_pdf': '[qcResDir]/[SeqID].insert.size.histogram.pdf'}
OUTPUT_KEYS = ['insert_size_metrics_txt', 'insert_size_hist_pdf']

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
    for k in ['BamDir', 'ReferenceFasta', 'SeqID', 'Threads', 'bind', 'gatk_jar', 'insert_size_hist_pdf', 'insert_size_metrics_txt', 'java_bin', 'lc_all', 'qcResDir', 'sif', 'xmx_mb']:
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
