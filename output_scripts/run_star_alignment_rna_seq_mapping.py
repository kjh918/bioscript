#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = star_alignment
# VERSION = 2.7.11a
# THREADS = 12
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
    parser.add_argument('--aligned_bam', required=False, default='[BamDir]/[SeqID].Aligned.sortedByCoord.out.bam', help='No description (Default: [BamDir]/[SeqID].Aligned.sortedByCoord.out.bam)')
    parser.add_argument('--aligned_transcriptome_bam', required=False, default='[BamDir]/[SeqID].Aligned.toTranscriptome.out.bam', help='No description (Default: [BamDir]/[SeqID].Aligned.toTranscriptome.out.bam)')
    parser.add_argument('--gene_counts', required=False, default='[BamDir]/[SeqID].ReadsPerGene.out.tab', help='No description (Default: [BamDir]/[SeqID].ReadsPerGene.out.tab)')
    parser.add_argument('--mapping_log', required=False, default='[BamDir]/[SeqID].Log.final.out', help='No description (Default: [BamDir]/[SeqID].Log.final.out)')
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
    aligned_bam = args.aligned_bam
    aligned_transcriptome_bam = args.aligned_transcriptome_bam
    gene_counts = args.gene_counts
    mapping_log = args.mapping_log
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
    if not aligned_bam:
    aligned_bam = f"{BamDir}/{SeqID}.Aligned.sortedByCoord.out.bam"
    if not aligned_transcriptome_bam:
    aligned_transcriptome_bam = f"{BamDir}/{SeqID}.Aligned.toTranscriptome.out.bam"
    if not gene_counts:
    gene_counts = f"{BamDir}/{SeqID}.ReadsPerGene.out.tab"
    if not mapping_log:
    mapping_log = f"{BamDir}/{SeqID}.Log.final.out"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {star_sif} STAR --runThreadN {Threads} --genomeDir {StarIndex} --readFilesIn {FastqDir}/{SeqID}{InputSuffix}_R1.fastq.gz {FastqDir}/{SeqID}{InputSuffix}_R2.fastq.gz --readFilesCommand zcat --quantMode {quantMode} --twopassMode {twopassMode} --chimSegmentMin {chimSegmentMin} --outFilterMismatchNmax {outFilterMismatchNmax} --outFileNamePrefix {BamDir}/{SeqID}. --outSAMtype {outSAMtype} --outSAMattributes {outSAMattributes} --outSAMunmapped {outSAMunmapped}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if aligned_transcriptome_bam:
        _tgt = os.path.dirname(aligned_transcriptome_bam) if os.path.splitext(aligned_transcriptome_bam)[1] else aligned_transcriptome_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if aligned_bam:
        _tgt = os.path.dirname(aligned_bam) if os.path.splitext(aligned_bam)[1] else aligned_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if FastqDir:
        _tgt = os.path.dirname(FastqDir) if os.path.splitext(FastqDir)[1] else FastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if mapping_log:
        _tgt = os.path.dirname(mapping_log) if os.path.splitext(mapping_log)[1] else mapping_log
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if gene_counts:
        _tgt = os.path.dirname(gene_counts) if os.path.splitext(gene_counts)[1] else gene_counts
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()