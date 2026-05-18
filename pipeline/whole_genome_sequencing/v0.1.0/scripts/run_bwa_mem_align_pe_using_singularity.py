#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bwa_mem
# VERSION = 0.7.17
# THREADS = 1
# PROFILE = align_pe_using_singularity

"""
Tool: bwa_mem (0.7.17)
Profile: align_pe_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bwa_mem Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming and SM tag (Default: )')
    parser.add_argument('--TrimFastqDir', required=True, default='', help='Directory containing FASTQ files (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Output directory for aligned BAM files (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--ReadGroupID', required=True, default='', help='Read Group ID (ID tag) (Default: )')
    parser.add_argument('--ReadGroupPlatform', required=True, default='', help='Sequencing platform (PL tag: e.g., ILLUMINA) (Default: )')
    parser.add_argument('--ReadGroupLibrary', required=True, default='', help='Library identifier (LB tag) (Default: )')
    parser.add_argument('--ReadGroupCenter', required=True, default='', help='Sequencing center name (CN tag) (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='trimmed', help='Suffix of input FASTQs (e.g., trimmed, raw, filtered) (Default: trimmed)')
    parser.add_argument('--OutputSuffix', required=False, default='primary', help='Suffix for output BAM (e.g., primary, initial) (Default: primary)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--bwa_bin', required=False, default='bwa', help='BWA binary path inside container (Default: bwa)')
    parser.add_argument('--samtools_bin', required=False, default='samtools', help='SAMtools binary path (host or inside if mounted) (Default: samtools)')
    parser.add_argument('--sif', required=False, default='/storage/images/bwa-0.7.17.sif', help='BWA singularity image (.sif) (Default: /storage/images/bwa-0.7.17.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='8', help='Number of threads for BWA and SAMtools (Default: 8)')
    parser.add_argument('--mark_short_split', required=False, default='-M', help='Mark shorter split hits as secondary (for Picard compatibility) (Default: -M)')
    parser.add_argument('--soft_clipping', required=False, default='-Y', help='Use soft clipping for supplementary alignments (Default: -Y)')
    parser.add_argument('--clipping_penalty', required=False, default='-L 50,50', help='Penalty for 5\'- and 3\'-end clipping (Default: -L 50,50)')
    parser.add_argument('--other_args', required=False, default='', help='Other additional BWA MEM arguments (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    TrimFastqDir = args.TrimFastqDir
    BamDir = args.BamDir
    ReferenceFasta = args.ReferenceFasta
    ReadGroupID = args.ReadGroupID
    ReadGroupPlatform = args.ReadGroupPlatform
    ReadGroupLibrary = args.ReadGroupLibrary
    ReadGroupCenter = args.ReadGroupCenter
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    bwa_bin = args.bwa_bin
    samtools_bin = args.samtools_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    mark_short_split = args.mark_short_split
    soft_clipping = args.soft_clipping
    clipping_penalty = args.clipping_penalty
    other_args = args.other_args

    # --- [Output Paths] ---
    primary_bam = f"{BamDir}/{SeqID}.{OutputSuffix}.bam"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {bwa_bin} mem  {mark_short_split} {soft_clipping} {clipping_penalty} {other_args}  -t {Threads} -R '@RG\tID:{ReadGroupID}\tPL:{ReadGroupPlatform}\tLB:{ReadGroupLibrary}\tSM:{SeqID}\tCN:{ReadGroupCenter}' {ReferenceFasta} {TrimFastqDir}/{SeqID}.{InputSuffix}_R1.fastq.gz {TrimFastqDir}/{SeqID}.{InputSuffix}_R2.fastq.gz | {samtools_bin} view -@ {Threads} -bS -o {primary_bam} -"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(TrimFastqDir) if '.' in os.path.basename(TrimFastqDir) else TrimFastqDir, exist_ok=True)
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()