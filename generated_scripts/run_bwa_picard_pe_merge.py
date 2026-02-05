#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = '[java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] FastqToSam --FASTQ [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz --FASTQ2 [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz --SAMPLE_NAME [SeqID] --OUTPUT [BamDir]/[SeqID].fastqtosam.bam --READ_GROUP_NAME [ReadGroupID] --PLATFORM [ReadGroupPlatform] --LIBRARY_NAME [ReadGroupLibrary] --SEQUENCING_CENTER [ReadGroupCenter] --TMP_DIR [TmpDir] && singularity exec -B [bind] [bwa_sif] bwa mem [bwa_args] -t [Threads] -R "@RG\\tID:[ReadGroupID]\\tPL:[ReadGroupPlatform]\\tLB:[ReadGroupLibrary]\\tSM:[SeqID]\\tCN:[ReadGroupCenter]" [ReferenceFasta] [TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz [TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz > [BamDir]/[SeqID].bwa.mem.sam && [java_bin] -XX:ParallelGCThreads=[gc_threads] -Xmx[xmx_mb]m -jar [picard_jar] MergeBamAlignment --UNMAPPED_BAM [BamDir]/[SeqID].fastqtosam.bam --ALIGNED_BAM [BamDir]/[SeqID].bwa.mem.sam --REFERENCE_SEQUENCE [ReferenceFasta] --OUTPUT [BamDir]/[SeqID].primary.bam [mba_args]'
REQUIRED_KEYS = ['SeqID', 'TrimFastqDir', 'BamDir', 'TmpDir', 'ReferenceFasta', 'ReadGroupID', 'ReadGroupPlatform', 'ReadGroupLibrary', 'ReadGroupCenter']
DEFAULTS = {'java_bin': 'java', 'picard_jar': '/storage/apps/bin/picard.jar', 'gc_threads': '14', 'xmx_mb': '16384', 'bind': '/storage,/data', 'bwa_sif': '/storage/images/bwa-0.7.17.sif', 'Threads': '8', 'bwa_args': '-M -Y -L 50,50', 'mba_args': '--CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF\n', 'unmapped_bam': '[BamDir]/[SeqID].fastqtosam.bam', 'aligned_sam': '[BamDir]/[SeqID].bwa.mem.sam', 'primary_bam': '[BamDir]/[SeqID].primary.bam', 'primary_bai': '[BamDir]/[SeqID].primary.bai'}
OUTPUT_KEYS = ['unmapped_bam', 'aligned_sam', 'primary_bam', 'primary_bai']

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
    for k in ['BamDir', 'ReadGroupCenter', 'ReadGroupID', 'ReadGroupLibrary', 'ReadGroupPlatform', 'ReferenceFasta', 'SeqID', 'Threads', 'TmpDir', 'TrimFastqDir', 'aligned_sam', 'bind', 'bwa_args', 'bwa_sif', 'gc_threads', 'java_bin', 'mba_args', 'picard_jar', 'primary_bai', 'primary_bam', 'unmapped_bam', 'xmx_mb']:
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
