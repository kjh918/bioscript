#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'mkdir -p [WorkDir] && ln -sf [NormalDir]/*.bed.gz [WorkDir]/ && ln -sf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/ &&\nfind [WorkDir] -maxdepth 1 -name "*.gz" | xargs -I {] bash -c \'N=$(basename "{]" | sed "s/\\.bed.*//"); [BinUnsorted] [BinFile] [BinCount] <(zcat "{]" | awk "{if(\\$1!~/^chr/) print \\"chr\\"\\$0; else print \\$0]") $N "{]"_mapped\' &&\nfind [WorkDir] -maxdepth 1 -name "*.bed" ! -name "*.gz" | xargs -I {] bash -c \'N=$(basename "{]" | sed "s/\\.bed.*//"); [BinUnsorted] [BinFile] [BinCount] <(awk "{if(\\$1!~/^chr/) print \\"chr\\"\\$0; else print \\$0]" "{]") $N "{]"_mapped\' &&\npaste [WorkDir]/*_mapped > [WorkDir]/data &&\n[ProcessR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Base] && [ReclustR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Reclust] && [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2'
REQUIRED_KEYS = ['SeqID', 'BinSize', 'ReadLength']
DEFAULTS = {'NormalSamples': 'Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz', 'R_Args_Base': 'status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1', 'R_Args_Reclust': 'status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0', 'BinMeth': 'variable_[BinSize]000_[ReadLength]_bwa', 'WorkDir': '[ResultBaseDir]/[SeqID]/[SeqID]_[BinSize]kb', 'BinUnsorted': '[GinkgoHomeDir]/scripts/binUnsorted', 'BinFile': '[GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth]', 'BinCount': '$(wc -l < [BinFile])', 'ProcessR': '[GinkgoHomeDir]/scripts/process.R', 'ReclustR': '[GinkgoHomeDir]/scripts/reclust.R', 'CNVCaller': '[GinkgoHomeDir]/scripts/CNVcaller'}
OUTPUT_KEYS = ['BinMeth', 'WorkDir', 'BinUnsorted', 'BinFile', 'BinCount', 'ProcessR', 'ReclustR', 'CNVCaller']

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
    for k in ['BamToBedDir', 'BinCount', 'BinFile', 'BinMeth', 'BinSize', 'BinUnsorted', 'CNVCaller', 'Genome', 'GinkgoHomeDir', 'NormalDir', 'NormalSamples', 'ProcessR', 'R_Args_Base', 'R_Args_Reclust', 'ReadLength', 'ReclustR', 'ResultBaseDir', 'SeqID', 'Threads', 'WorkDir']:
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
