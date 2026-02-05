#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'mkdir -p [AnalysisRunDir] && ln -Tsf [NGS_DataBaseDir]/[SeqID].analysisReady.bam [AnalysisRunDir]/[SampleID].bam && ln -Tsf [NGS_DataBaseDir]/[SeqID].analysisReady.bam.bai [AnalysisRunDir]/[SampleID].bam.bai && Rscript [Rscript_path] \\ --SeqID [SampleID] \\ --AnalysisRunDir [AnalysisRunDir] \\ --BinSize [BinSize] \\ --Ploidy [SamplePloidy] \\  --GenomeVersion [GenomeVersion]'
REQUIRED_KEYS = ['SeqID', 'SampleID', 'SamplePloidy', 'NGS_DataBaseDir', 'AnalysisRunDir']
DEFAULTS = {'Threads': '1', 'BinSize': '200kb', 'GenomeVersion': 'hg38', 'Rscript_path': '/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R', 'SamplePloidy': '2', 'AnalysisRunDir': '[ResultBaseDir]/[SampleID]'}
OUTPUT_KEYS = ['AnalysisRunDir']

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
    for k in ['AnalysisRunDir', 'BinSize', 'GenomeVersion', 'NGS_DataBaseDir', 'Rscript_path', 'SampleID', 'SamplePloidy', 'SeqID', 'Threads']:
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
