#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'mkdir -p [TmpDir] && \n[GatkPath] --java-options "-XX:ParallelGCThreads=[Threads] -Xmx[Memory] -Djava.io.tmpdir=[TmpDir]"  HaplotypeCaller  -R [ReferenceFasta] -I [InputBam] -L [Chromosome]  -ploidy [Ploidy] -stand-call-conf 30 --dbsnp [DbsnpVcf]  -O [OutGvcf] -ERC GVCF && \n[GatkPath] --java-options "-XX:ParallelGCThreads=[Threads] -Xmx[Memory] -Djava.io.tmpdir=[TmpDir]"  GenotypeGVCFs  --include-non-variant-sites [IncludeNonVariant]  -R [ReferenceFasta] -V [OutGvcf] -O [OutVcf]'
REQUIRED_KEYS = ['SeqID', 'Chromosome', 'InputBam', 'ResultDir']
DEFAULTS = {'GatkPath': '/storage/apps/gatk-4.4.0.0/gatk', 'Threads': '4', 'Memory': '32g', 'Ploidy': '2', 'TmpDir': '[ResultDir]/tmp', 'IncludeNonVariant': 'false', 'OutVcf': '[ResultDir]/[SeqID].[Chromosome].vcf.gz', 'OutGvcf': '[ResultDir]/[SeqID].[Chromosome].gvcf.gz'}
OUTPUT_KEYS = ['OutVcf', 'OutGvcf']

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
    for k in ['Chromosome', 'DbsnpVcf', 'GatkPath', 'IncludeNonVariant', 'InputBam', 'Memory', 'OutGvcf', 'OutVcf', 'Ploidy', 'ReferenceFasta', 'ResultDir', 'SeqID', 'Threads', 'TmpDir']:
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
