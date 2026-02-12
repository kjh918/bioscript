#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = '[BcftoolsPath] mpileup -f [ReferenceFasta] -T [SitesVcfGz] -q [MinMQ] -Q [MinBQ] -a FORMAT/AD,FORMAT/DP -Ou [InputBamDir]/[SeqID].analysisReady.bam | [BcftoolsPath] call -Am -Oz -o [RawVcf] && [BcftoolsPath] index -f [RawVcf] && [BcftoolsPath] annotate -a [PopAfAnnotVcf] -c [AnnotationQuery] -h [PopAfHeaderHdr] -Oz -o [AnnVcf] [RawVcf] && [BcftoolsPath] index -f [AnnVcf]'
REQUIRED_KEYS = ['SeqID', 'Chromosome', 'InputBamDir', 'ReferenceFasta', 'SitesVcfGz', 'PopAfAnnotVcf', 'PopAfHeaderHdr', 'ResultDir']
DEFAULTS = {'BcftoolsPath': '/storage/home/jhkim/Apps/bcftools/bcftools', 'BgzipPath': 'bgzip', 'Threads': '4', 'MinBQ': '20', 'MinMQ': '30', 'TmpDir': '[ResultDir]/tmp', 'AnnotationQuery': 'CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN', 'VcfQuery': '%CHROM|%POS|%REF|%ALT|%INFO/KOVA_AF|%INFO/KOVA_AN|%INFO/GNOMAD_AF|%INFO/GNOMAD_AN|[%GT]|[%DP]|[%AD]\\n', 'RawVcf': '[ResultDir]/[SeqID].[Chromosome].sites.raw.vcf.gz', 'AnnVcf': '[ResultDir]/[SeqID].[Chromosome].sites.af.vcf.gz'}
OUTPUT_KEYS = ['RawVcf', 'AnnVcf']

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
    for k in ['AnnVcf', 'AnnotationQuery', 'BcftoolsPath', 'BgzipPath', 'Chromosome', 'InputBamDir', 'MinBQ', 'MinMQ', 'PopAfAnnotVcf', 'PopAfHeaderHdr', 'RawVcf', 'ReferenceFasta', 'ResultDir', 'SeqID', 'SitesVcfGz', 'Threads', 'TmpDir', 'VcfQuery']:
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
