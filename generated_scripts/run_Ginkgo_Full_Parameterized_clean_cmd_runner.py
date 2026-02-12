#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'for n in [NormalSamples]; do ln -Tsf [NormalDir]/[n] [WorkDir]/[n]; done && ln -Tsf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/[SeqID].bed.gz && ls [WorkDir] | grep ".bed.gz$" > [WorkDir]/list &&\nwhile read -r f; do\n  echo ">> Binning: [f]" &&\n  SampleName="[f%.bed.gz]" &&\n  SampleName="[SampleName%.bed]" &&\n  \n  # [Normalization] chr 접두사 보정\n  if [[ "[f]" =~ \\.gz$ ]]; then\n    zcat [WorkDir]/[f] | awk \'{if($1!~/^chr/) print "chr"$0; else print $0]\' | gzip > [WorkDir]/[f]_tmp.gz &&\n    [BinUnsorted] [BinFile] [BinCount] <(zcat [WorkDir]/[f]_tmp.gz) [SampleName] [WorkDir]/[f]_mapped;\n  else\n    awk \'{if($1!~/^chr/) print "chr"$0; else print $0]\' [WorkDir]/[f] > [WorkDir]/[f]_tmp &&\n    [BinUnsorted] [BinFile] [BinCount] [WorkDir]/[f]_tmp [SampleName] [WorkDir]/[f]_mapped;\n  fi;\ndone < [WorkDir]/list &&\n# 데이터 통합 및 분석 실행 paste [WorkDir]/*_mapped > [WorkDir]/data && cd [GinkgoHomeDir] && [ProcessR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Base] && [ReclustR] [GinkgoHomeDir]/genomes/[Genome]/original [WorkDir] [R_Args_Reclust] && [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2 &&\n# 결과 정리 mkdir -p [FinalDir] && mv [WorkDir]/* [FinalDir]/ && rm -rf [WorkDir]'
REQUIRED_KEYS = ['SeqID', 'BinSize', 'ReadLength']
DEFAULTS = {'NormalSamples': 'Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz', 'R_Args_Base': 'status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1', 'R_Args_Reclust': 'status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0', 'BinMeth': 'variable_[BinSize]000_[ReadLength]_bwa', 'WorkDir': '[GinkgoHomeDir]/uploads/[SeqID]_[BinSize]kb', 'FinalDir': '[ResultBaseDir]/[SeqID]/Binsize_[BinSize]kb', 'BinUnsorted': '[GinkgoHomeDir]/scripts/binUnsorted', 'BinFile': '[GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth]', 'BinCount': '$(wc -l < [BinFile])', 'ProcessR': '[GinkgoHomeDir]/scripts/process.R', 'ReclustR': '[GinkgoHomeDir]/scripts/reclust.R', 'CNVCaller': '[GinkgoHomeDir]/scripts/CNVcaller'}
OUTPUT_KEYS = ['BinMeth', 'WorkDir', 'FinalDir', 'BinUnsorted', 'BinFile', 'BinCount', 'ProcessR', 'ReclustR', 'CNVCaller']

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
    for k in ['BamToBedDir', 'BinCount', 'BinFile', 'BinMeth', 'BinSize', 'BinUnsorted', 'CNVCaller', 'FinalDir', 'Genome', 'GinkgoHomeDir', 'NormalDir', 'NormalSamples', 'Normalization', 'ProcessR', 'R_Args_Base', 'R_Args_Reclust', 'ReadLength', 'ReclustR', 'ResultBaseDir', 'SampleName', 'SeqID', 'WorkDir', 'f', 'n']:
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
