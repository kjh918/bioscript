#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'samtools view -b -h -q [min_mapq] -f [include_flag] -F [exclude_flag1] -F [exclude_flag2] -e "[expr_sclen]" -e "[expr_nm]" --threads [Threads] [BamDir]/[SeqID].recal.bam > [BamDir]/[SeqID].recal.filtered.bam && samtools index --threads [Threads] [BamDir]/[SeqID].recal.filtered.bam && ln -Tsf  [BamDir]/[SeqID].recal.filtered.bam [BamDir]/[SeqID].analysisReady.bam && ln -Tsf  [BamDir]/[SeqID].recal.filtered.bam.bai [BamDir]/[SeqID].analysisReady.bam.bai'
REQUIRED_KEYS = ['SeqID', 'BamDir']
DEFAULTS = {'Threads': '8', 'min_mapq': '20', 'include_flag': '0x2', 'exclude_flag1': '0x100', 'exclude_flag2': '0x4', 'expr_sclen': 'sclen < 20', 'expr_nm': '[NM] < 12', 'filtered_bam': '[BamDir]/[SeqID].recal.filtered.bam', 'filtered_bai': '[BamDir]/[SeqID].recal.filtered.bam.bai', 'analysis_ready_bam': '[BamDir]/[SeqID].analysisReady.bam', 'analysis_ready_bai': '[BamDir]/[SeqID].analysisReady.bam.bai'}
OUTPUT_KEYS = ['filtered_bam', 'filtered_bai', 'analysis_ready_bam', 'analysis_ready_bai']

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
    for k in ['BamDir', 'SeqID', 'Threads', 'analysis_ready_bai', 'analysis_ready_bam', 'exclude_flag1', 'exclude_flag2', 'expr_nm', 'expr_sclen', 'filtered_bai', 'filtered_bam', 'include_flag', 'min_mapq']:
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
