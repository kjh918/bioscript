#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'singularity exec -B [bind] [sif] /opt/mosdepth --threads [Threads] [mosdepth_args] --by [bin] --mapq [mapq] [qcResDir]/[SeqID] [BamDir]/[SeqID].analysisReady.bam'
REQUIRED_KEYS = ['SeqID', 'BamDir', 'qcResDir']
DEFAULTS = {'bind': '/storage,/data', 'sif': '/storage/images/mosdepth-0.3.6.sif', 'Threads': '8', 'bin': '100000', 'mapq': '20', 'mosdepth_args': '--no-per-base --fast-mode', 'prefix': '[qcResDir]/[SeqID]', 'regions_bed_gz': '[qcResDir]/[SeqID].regions.bed.gz', 'regions_bed_gz_csi': '[qcResDir]/[SeqID].regions.bed.gz.csi', 'summary_txt': '[qcResDir]/[SeqID].mosdepth.summary.txt', 'global_dist_txt': '[qcResDir]/[SeqID].global.dist.txt'}
OUTPUT_KEYS = ['prefix', 'regions_bed_gz', 'regions_bed_gz_csi', 'summary_txt', 'global_dist_txt']

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
    for k in ['BamDir', 'SeqID', 'Threads', 'bin', 'bind', 'global_dist_txt', 'mapq', 'mosdepth_args', 'prefix', 'qcResDir', 'regions_bed_gz', 'regions_bed_gz_csi', 'sif', 'summary_txt']:
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
