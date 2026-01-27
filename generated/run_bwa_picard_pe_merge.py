#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict

TOKEN_PAT = re.compile(r"\[([A-Za-z_][A-Za-z0-9_]*)\]")

CMD_LINE = '[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] FastqToSam --FASTQ [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz --FASTQ2 [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz --SAMPLE_NAME [SeqID] --OUTPUT [BamDir]/[SeqID].fastqtosam.bam --READ_GROUP_NAME [ReadGroupID] --PLATFORM [ReadGroupPlatform] --LIBRARY_NAME [ReadGroupLibrary] --SEQUENCING_CENTER [ReadGroupCenter] --TMP_DIR [TmpDir] && singularity exec -B [bind] [bwa_sif] bwa mem [bwa_args] -t [Threads] -R "@RG\\tID:[ReadGroupID]\\tPL:[ReadGroupPlatform]\\tLB:[ReadGroupLibrary]\\tSM:[SeqID]\\tCN:[ReadGroupCenter]" [ReferenceFasta] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz > [BamDir]/[SeqID].bwa.mem.sam && [java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] MergeBamAlignment --UNMAPPED_BAM [BamDir]/[SeqID].fastqtosam.bam --ALIGNED_BAM [BamDir]/[SeqID].bwa.mem.sam --REFERENCE_SEQUENCE [ReferenceFasta] --OUTPUT [BamDir]/[SeqID].primary.bam [mba_args]'
ALL_KEYS = ['SeqID', 'TrimFastqDir', 'BamDir', 'TmpDir', 'ReferenceFasta', 'ReadGroupID', 'ReadGroupPlatform', 'ReadGroupLibrary', 'ReadGroupCenter', 'java_bin', 'picard_jar', 'gc_threads', 'xmx_mb', 'bind', 'bwa_sif', 'Threads', 'bwa_args', 'mba_args', 'unmapped_bam', 'aligned_sam', 'primary_bam', 'primary_bai']
REQUIRED = {'SeqID': True, 'TrimFastqDir': True, 'BamDir': True, 'TmpDir': True, 'ReferenceFasta': True, 'ReadGroupID': True, 'ReadGroupPlatform': True, 'ReadGroupLibrary': True, 'ReadGroupCenter': True, 'unmapped_bam': False, 'aligned_sam': False, 'primary_bam': False, 'primary_bai': False}
DEFAULTS = {'java_bin': 'java', 'picard_jar': '/storage/apps/bin/picard.jar', 'gc_threads': '14', 'xmx_mb': '16384', 'bind': '/storage,/data', 'bwa_sif': '/storage/images/bwa-0.7.17.sif', 'Threads': '8', 'bwa_args': '-M -Y -L 50,50', 'mba_args': '--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF\n', 'unmapped_bam': '[BamDir]/[SeqID].fastqtosam.bam', 'aligned_sam': '[BamDir]/[SeqID].bwa.mem.sam', 'primary_bam': '[BamDir]/[SeqID].primary.bam', 'primary_bai': '[BamDir]/[SeqID].primary.bai'}
OUTPUTS_META = {'unmapped_bam': {'default': '[BamDir]/[SeqID].fastqtosam.bam'}, 'aligned_sam': {'default': '[BamDir]/[SeqID].bwa.mem.sam'}, 'primary_bam': {'default': '[BamDir]/[SeqID].primary.bam'}, 'primary_bai': {'default': '[BamDir]/[SeqID].primary.bai'}}

def render(s: str, ctx: Dict[str, str]) -> str:
    def repl(m: re.Match) -> str:
        key = m.group(1)
        return (ctx.get(key) or "").strip()
    out = TOKEN_PAT.sub(repl, s)
    out = re.sub(r"\s+", " ", out).strip()
    return out

def finalize_outputs(ctx: Dict[str, str]) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    for name, meta in OUTPUTS_META.items():
        meta = meta or {}
        cur = (ctx.get(name) or "").strip()
        if cur:
            outputs[name] = cur
            continue

        dv = meta.get("default")
        if isinstance(dv, str) and dv.strip():
            val = render(dv, ctx)
            outputs[name] = val
            ctx[name] = val
            continue

        if bool(meta.get("required") is True):
            raise SystemExit(f"Missing required output --{name} (no default/template inferred)")
    return outputs

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Runner: bwa_picard bwa_0.7.17-picard_3.1.0 (pe_merge)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    for k in ALL_KEYS:
        if REQUIRED.get(k, False) and k not in OUTPUTS_META:
            ap.add_argument(f"--{k}", required=True, default="")
        else:
            ap.add_argument(f"--{k}", required=False, default=DEFAULTS.get(k, ""))

    ap.add_argument("--cwd", default=".")
    ap.add_argument("--print-cmd", action="store_true")
    ap.add_argument("--print-outputs", action="store_true")
    ap.add_argument("--emit-outputs", default="", help="write outputs json to file path")

    args = ap.parse_args()
    ctx = {k: ("" if v is None else str(v)) for k, v in vars(args).items()
            if k not in ("cwd","print_cmd","print_outputs","emit_outputs")}

    for k, dv in DEFAULTS.items():
        if not ctx.get(k, "").strip():
            ctx[k] = str(dv)

    outputs = finalize_outputs(ctx)

    cmd = render(CMD_LINE, ctx)

    if args.print_cmd:
        print(cmd)
        return

    proc = subprocess.run(shlex.split(cmd), cwd=str(Path(args.cwd)))
    rc = proc.returncode

    if args.print_outputs:
        print(json.dumps(outputs, ensure_ascii=False))

    if args.emit_outputs:
        Path(args.emit_outputs).write_text(json.dumps(outputs, ensure_ascii=False, indent=2), encoding="utf-8")

    raise SystemExit(rc)

if __name__ == "__main__":
    main()
