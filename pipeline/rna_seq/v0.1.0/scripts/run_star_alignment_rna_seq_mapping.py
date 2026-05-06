#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = star_alignment
# VERSION = 2.7.11a
# THREADS = 1
# PROFILE = rna_seq_mapping

"""
Tool: star_alignment (2.7.11a)
Profile: rna_seq_mapping
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="star_alignment Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--FastqDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--StarIndex', required=True, default='', help='Path to STAR index directory (Default: )')
    parser.add_argument('--Threads', required=False, default='12', help='No description (Default: 12)')
    parser.add_argument('--InputSuffix', required=False, default='.trimmed', help='No description (Default: .trimmed)')
    parser.add_argument('--outSAMtype', required=False, default='BAM SortedByCoordinate', help='No description (Default: BAM SortedByCoordinate)')
    parser.add_argument('--quantMode', required=False, default='TranscriptomeSAM', help='No description (Default: TranscriptomeSAM)')
    parser.add_argument('--outSAMattributes', required=False, default='NH HI AS nM XS NM', help='No description (Default: NH HI AS nM XS NM)')
    parser.add_argument('--chimSegmentMin', required=False, default='10', help='No description (Default: 10)')
    parser.add_argument('--twopassMode', required=False, default='Basic', help='No description (Default: Basic)')
    parser.add_argument('--outFilterMismatchNmax', required=False, default='10', help='No description (Default: 10)')
    parser.add_argument('--outSAMunmapped', required=False, default='Within', help='No description (Default: Within)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--star_sif', required=False, default='/storage/images/star-2.7.11.sif', help='No description (Default: /storage/images/star-2.7.11.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    FastqDir = args.FastqDir
    BamDir = args.BamDir
    StarIndex = args.StarIndex
    Threads = args.Threads
    InputSuffix = args.InputSuffix
    outSAMtype = args.outSAMtype
    quantMode = args.quantMode
    outSAMattributes = args.outSAMattributes
    chimSegmentMin = args.chimSegmentMin
    twopassMode = args.twopassMode
    outFilterMismatchNmax = args.outFilterMismatchNmax
    outSAMunmapped = args.outSAMunmapped
    singularity_bin = args.singularity_bin
    star_sif = args.star_sif
    bind = args.bind

    # --- [Output Paths] ---
    aligned_bam = f"{BamDir}/{SeqID}.Aligned.sortedByCoord.out.bam"
    gene_counts = f"{BamDir}/{SeqID}.ReadsPerGene.out.tab"
    mapping_log = f"{BamDir}/{SeqID}.Log.final.out"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {star_sif} STAR --runThreadN {Threads} --genomeDir {StarIndex} --readFilesIn {FastqDir}/{SeqID}{InputSuffix}_R1.fastq.gz {FastqDir}/{SeqID}{InputSuffix}_R2.fastq.gz --readFilesCommand zcat --quantMode {quantMode} --twopassMode {twopassMode} --chimSegmentMin {chimSegmentMin} --outFilterMismatchNmax {outFilterMismatchNmax} --outFileNamePrefix {BamDir}/{SeqID}. --outSAMtype {outSAMtype} --outSAMattributes {outSAMattributes} --outSAMunmapped {outSAMunmapped}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(FastqDir) if '.' in os.path.basename(FastqDir) else FastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()