#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'echo "==== [Integrated Ginkgo] Starting [RunID] ====" && mkdir -p [WorkDir] &&\n# 1. 데이터 준비 및 링크 for normal in [NormalSamples]; do ln -Tsf [NormalDir]/[normal] [WorkDir]/[normal]; done && ln -Tsf [BamToBedDir]/[SeqID].bed.gz [WorkDir]/[SeqID].bed.gz && ls [WorkDir] | grep ".bed.gz$" > [WorkDir]/list &&\n# 2. 전처리 (analyze.sh 로직 재현: chr 접두사 체크 및 Binning) while read file; do \\\n  echo ">> Binning: [file]"; \\\n  if [[ "[file]" =~ \\.gz$ ]]; then \\\n    zcat [WorkDir]/[file] | awk \'{if($1!~/^chr/) print "chr"$0; else print $0]\' | gzip > [WorkDir]/[file]_tmp.gz; \\\n    [BinUnsorted] [GenomeDir]/[BinMeth] $(wc -l < [GenomeDir]/[BinMeth]) <(zcat [WorkDir]/[file]_tmp.gz) [file%.bed.gz] [WorkDir]/[file]_mapped; \\\n  else \\\n    awk \'{if($1!~/^chr/) print "chr"$0; else print $0]\' [WorkDir]/[file] > [WorkDir]/[file]_tmp; \\\n    [BinUnsorted] [GenomeDir]/[BinMeth] $(wc -l < [GenomeDir]/[BinMeth]) [WorkDir]/[file]_tmp [file%.bed] [WorkDir]/[file]_mapped; \\\n  fi; \\\ndone < [WorkDir]/list &&\n# 3. 데이터 통합 (Merge binned reads) paste [WorkDir]/*_mapped > [WorkDir]/data &&\n# 4. R 분석 실행 (Ginkgo Core) cd [GinkgoHomeDir] && \\ [ProcessR] [GenomeDir] [WorkDir] status.xml data 2 [BinMeth] [ClustMeth] [DistMeth] 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1 && \\ [ReclustR] [GenomeDir] [WorkDir] status.xml [BinMeth] [ClustMeth] [DistMeth] 0 ploidyDummy.txt 0 &&\n# 5. CNV Calling (analyze.sh 후반부 로직) [CNVCaller] [WorkDir]/SegCopy [WorkDir]/CNV1 [WorkDir]/CNV2 &&\n# 6. 결과 이동 및 청소 mkdir -p [FinalDir] && \\ mv [WorkDir]/* [FinalDir]/ && \\ rm -rf [WorkDir] && \\ echo "==== [Integrated Ginkgo] Completed [RunID] ===="'
REQUIRED_KEYS = ['SeqID', 'BinSize', 'ReadLength']
DEFAULTS = {'NormalSamples': 'Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz', 'BinMeth': 'variable_[BinSize]000_[ReadLength]_bwa', 'RunID': '[SeqID]_[BinSize]kb', 'WorkDir': '[GinkgoHomeDir]/uploads/[SeqID]_[BinSize]kb', 'FinalDir': '[ResultBaseDir]/[SeqID]/Binsize_[BinSize]kb', 'BinUnsorted': '[GinkgoHomeDir]/scripts/binUnsorted', 'ProcessR': '[GinkgoHomeDir]/scripts/process.R', 'ReclustR': '[GinkgoHomeDir]/scripts/reclust.R', 'CNVCaller': '[GinkgoHomeDir]/scripts/CNVcaller', 'GenomeDir': '[GinkgoHomeDir]/genomes/[Genome]/original'}
OUTPUT_KEYS = ['RunID', 'WorkDir', 'FinalDir', 'BinUnsorted', 'ProcessR', 'ReclustR', 'CNVCaller', 'GenomeDir']

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
    for k in ['BamToBedDir', 'BinMeth', 'BinSize', 'BinUnsorted', 'CNVCaller', 'ClustMeth', 'DistMeth', 'FinalDir', 'Genome', 'GenomeDir', 'GinkgoHomeDir', 'NormalDir', 'NormalSamples', 'ProcessR', 'ReadLength', 'ReclustR', 'ResultBaseDir', 'RunID', 'SeqID', 'WorkDir', 'file', 'normal']:
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
