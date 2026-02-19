#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = qc_insert_size_using_singularity

"""
Tool: gatk4 (4.4.0.0)
Profile: qc_insert_size_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the BAM file to be analyzed (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Output directory for QC metrics and PDF histogram (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of the input BAM file (e.g., analysisReady, recal, primary) (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--java_bin', required=False, default='java', help='Path to java executable inside or outside container (Default: java)')
    parser.add_argument('--gatk_jar', required=False, default='/gatk/gatk-package-4.4.0.0-local.jar', help='Path to GATK local jar file inside container (Default: /gatk/gatk-package-4.4.0.0-local.jar)')
    parser.add_argument('--sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='Path to GATK singularity image (.sif) (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='14', help='Number of threads for ParallelGC (Default: 14)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Maximum Java heap memory (MB) (Default: 16384)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    java_bin = args.java_bin
    gatk_jar = args.gatk_jar
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb

    # --- [Output Paths] ---
    insert_size_metrics_txt = f"{qcResDir}/{SeqID}.insert_size.metrics.txt"
    insert_size_hist_pdf = f"{qcResDir}/{SeqID}.insert_size.histogram.pdf"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {java_bin} -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m -jar {gatk_jar} CollectInsertSizeMetrics --INPUT {BamDir}/{SeqID}.bam --OUTPUT {insert_size_metrics_txt} --Histogram_FILE {insert_size_hist_pdf} --REFERENCE_SEQUENCE {ReferenceFasta}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()