#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'singularity exec -B [bind] [sif] fastqc [fastqc_args] --threads [threads] --outdir [qcResDir] [RawFastqDir]/[SeqID]_R1.fastq.gz [RawFastqDir]/[SeqID]_R2.fastq.gz'
REQUIRED_KEYS = ['SeqID', 'RawFastqDir', 'qcResDir']
DEFAULTS = {'bind': '/storage,/data', 'sif': '/storage/images/fastqc-0.12.1.sif', 'threads': '8', 'fastqc_args': '--extract', 'r1_html': '[qcResDir]/[SeqID]_R1_fastqc.html', 'r2_html': '[qcResDir]/[SeqID]_R2_fastqc.html', 'r1_zip': '[qcResDir]/[SeqID]_R1_fastqc.zip', 'r2_zip': '[qcResDir]/[SeqID]_R2_fastqc.zip', 'r1_dir': '[qcResDir]/[SeqID]_R1_fastqc', 'r2_dir': '[qcResDir]/[SeqID]_R2_fastqc'}
OUTPUT_KEYS = ['r1_html', 'r2_html', 'r1_zip', 'r2_zip', 'r1_dir', 'r2_dir']

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
    for k in ['RawFastqDir', 'SeqID', 'bind', 'fastqc_args', 'qcResDir', 'r1_dir', 'r1_html', 'r1_zip', 'r2_dir', 'r2_html', 'r2_zip', 'sif', 'threads']:
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
