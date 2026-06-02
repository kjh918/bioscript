#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bwa_picard
# VERSION = bwa_0.7.17-picard_3.1.0
# THREADS = 8
# PROFILE = pe_align_merge

"""
Tool: bwa_picard (bwa_0.7.17-picard_3.1.0)
Profile: pe_align_merge
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bwa_picard Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier for sample tracking (Default: )')
    parser.add_argument('--TrimFastqDir', required=True, default='', help='Directory containing input FASTQ files (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory to store intermediate and final BAM files (Default: )')
    parser.add_argument('--TmpDir', required=True, default='', help='Temporary directory for Picard operations (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--ReadGroupID', required=True, default='', help='Read Group ID (ID tag) (Default: )')
    parser.add_argument('--ReadGroupPlatform', required=True, default='', help='Sequencing platform (PL tag) (Default: )')
    parser.add_argument('--ReadGroupLibrary', required=True, default='', help='Library identifier (LB tag) (Default: )')
    parser.add_argument('--ReadGroupCenter', required=True, default='', help='Sequencing center (CN tag) (Default: )')
    parser.add_argument('--unmapped_bam', required=False, default='[BamDir]/[SeqID][InputSuffix].fastqtosam.bam', help='Intermediate unmapped BAM containing original metadata (Default: [BamDir]/[SeqID][InputSuffix].fastqtosam.bam)')
    parser.add_argument('--aligned_sam', required=False, default='[BamDir]/[SeqID][InputSuffix].bwa.mem.sam', help='Raw BWA alignment output without metadata (Default: [BamDir]/[SeqID][InputSuffix].bwa.mem.sam)')
    parser.add_argument('--primary_bam', required=False, default='[BamDir]/[SeqID][OutputSuffix].bam', help='Final merged and coordinate-sorted BAM file (Default: [BamDir]/[SeqID][OutputSuffix].bam)')
    parser.add_argument('--primary_bai', required=False, default='[BamDir]/[SeqID][OutputSuffix].bai', help='Index for the final merged BAM (Default: [BamDir]/[SeqID][OutputSuffix].bai)')
    parser.add_argument('--InputSuffix', required=False, default='.trimmed', help='Suffix of input FASTQs (e.g., trimmed, raw) (Default: .trimmed)')
    parser.add_argument('--OutputSuffix', required=False, default='.primary', help='Suffix for the final merged result (Default: .primary)')
    parser.add_argument('--java_bin', required=False, default='java', help='Path to Java executable (Default: java)')
    parser.add_argument('--picard_jar', required=False, default='/storage/apps/bin/picard.jar', help='Path to Picard.jar (Default: /storage/apps/bin/picard.jar)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to Singularity (Default: singularity)')
    parser.add_argument('--bwa_sif', required=False, default='/storage/images/bwa-0.7.17.sif', help='BWA Singularity image (Default: /storage/images/bwa-0.7.17.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for Singularity (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Max Java heap memory (MB) (Default: 16384)')
    parser.add_argument('--Threads', required=False, default='8', help='Threads for BWA and Java ParallelGC (Default: 8)')
    parser.add_argument('--bwa_bin', required=False, default='bwa', help='No description (Default: bwa)')
    parser.add_argument('--mark_short_split', required=False, default='-M', help='No description (Default: -M)')
    parser.add_argument('--soft_clipping', required=False, default='', help='No description (Default: )')
    parser.add_argument('--clipping_penalty', required=False, default='', help='No description (Default: )')
    parser.add_argument('--other_args', required=False, default='', help='No description (Default: )')
    parser.add_argument('--mba_strategy', required=False, default='MostDistant', help='No description (Default: MostDistant)')
    parser.add_argument('--mba_attributes', required=False, default='XS', help='No description (Default: XS)')
    parser.add_argument('--mba_orientations', required=False, default='FR --EXPECTED_ORIENTATIONS RF', help='No description (Default: FR --EXPECTED_ORIENTATIONS RF)')
    parser.add_argument('--mba_other_flags', required=False, default='--INCLUDE_SECONDARY_ALIGNMENTS true --ADD_MATE_CIGAR true --CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true', help='No description (Default: --INCLUDE_SECONDARY_ALIGNMENTS true --ADD_MATE_CIGAR true --CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    TrimFastqDir = args.TrimFastqDir
    BamDir = args.BamDir
    TmpDir = args.TmpDir
    ReferenceFasta = args.ReferenceFasta
    ReadGroupID = args.ReadGroupID
    ReadGroupPlatform = args.ReadGroupPlatform
    ReadGroupLibrary = args.ReadGroupLibrary
    ReadGroupCenter = args.ReadGroupCenter
    unmapped_bam = args.unmapped_bam
    aligned_sam = args.aligned_sam
    primary_bam = args.primary_bam
    primary_bai = args.primary_bai
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    java_bin = args.java_bin
    picard_jar = args.picard_jar
    singularity_bin = args.singularity_bin
    bwa_sif = args.bwa_sif
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads
    bwa_bin = args.bwa_bin
    mark_short_split = args.mark_short_split
    soft_clipping = args.soft_clipping
    clipping_penalty = args.clipping_penalty
    other_args = args.other_args
    mba_strategy = args.mba_strategy
    mba_attributes = args.mba_attributes
    mba_orientations = args.mba_orientations
    mba_other_flags = args.mba_other_flags

    # --- [Output Paths] ---
    if not unmapped_bam:
    unmapped_bam = f"{BamDir}/{SeqID}{InputSuffix}.fastqtosam.bam"
    if not aligned_sam:
    aligned_sam = f"{BamDir}/{SeqID}{InputSuffix}.bwa.mem.sam"
    if not primary_bam:
    primary_bam = f"{BamDir}/{SeqID}{OutputSuffix}.bam"
    if not primary_bai:
    primary_bai = f"{BamDir}/{SeqID}{OutputSuffix}.bai"

    # --- [Command Execution] ---
    cmd = f"{java_bin} -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m -jar {picard_jar} FastqToSam --FASTQ {TrimFastqDir}/{SeqID}{InputSuffix}_R1.fastq.gz --FASTQ2 {TrimFastqDir}/{SeqID}{InputSuffix}_R2.fastq.gz --SAMPLE_NAME {SeqID} --OUTPUT {unmapped_bam} --READ_GROUP_NAME {ReadGroupID} --PLATFORM {ReadGroupPlatform} --LIBRARY_NAME {ReadGroupLibrary} --SEQUENCING_CENTER {ReadGroupCenter} --TMP_DIR {TmpDir} && {singularity_bin} exec -B {bind} {bwa_sif} {bwa_bin} mem  {mark_short_split} {soft_clipping} {clipping_penalty} {other_args}  -t {Threads} -R '@RG\\tID:{ReadGroupID}\\tPL:{ReadGroupPlatform}\\tLB:{ReadGroupLibrary}\\tSM:{SeqID}\\tCN:{ReadGroupCenter}' {ReferenceFasta} {TrimFastqDir}/{SeqID}{InputSuffix}_R1.fastq.gz {TrimFastqDir}/{SeqID}{InputSuffix}_R2.fastq.gz > {aligned_sam} && {java_bin} -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m -jar {picard_jar} MergeBamAlignment --UNMAPPED_BAM {unmapped_bam} --ALIGNED_BAM {aligned_sam} --REFERENCE_SEQUENCE {ReferenceFasta} --OUTPUT {primary_bam} --PRIMARY_ALIGNMENT_STRATEGY {mba_strategy} --ATTRIBUTES_TO_RETAIN {mba_attributes} --EXPECTED_ORIENTATIONS {mba_orientations}  {mba_other_flags}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if primary_bai:
        _tgt = os.path.dirname(primary_bai) if os.path.splitext(primary_bai)[1] else primary_bai
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if primary_bam:
        _tgt = os.path.dirname(primary_bam) if os.path.splitext(primary_bam)[1] else primary_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if TmpDir:
        _tgt = os.path.dirname(TmpDir) if os.path.splitext(TmpDir)[1] else TmpDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if aligned_sam:
        _tgt = os.path.dirname(aligned_sam) if os.path.splitext(aligned_sam)[1] else aligned_sam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if TrimFastqDir:
        _tgt = os.path.dirname(TrimFastqDir) if os.path.splitext(TrimFastqDir)[1] else TrimFastqDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if unmapped_bam:
        _tgt = os.path.dirname(unmapped_bam) if os.path.splitext(unmapped_bam)[1] else unmapped_bam
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()