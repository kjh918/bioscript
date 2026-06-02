#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 14
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
    parser.add_argument('--insert_size_metrics_txt', required=False, default='[qcResDir]/[SeqID].insert_size.metrics.txt', help='Text file containing detailed insert size statistics (Default: [qcResDir]/[SeqID].insert_size.metrics.txt)')
    parser.add_argument('--insert_size_hist_pdf', required=False, default='[qcResDir]/[SeqID].insert_size.histogram.pdf', help='PDF file containing the insert size distribution histogram (Default: [qcResDir]/[SeqID].insert_size.histogram.pdf)')
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
    insert_size_metrics_txt = args.insert_size_metrics_txt
    insert_size_hist_pdf = args.insert_size_hist_pdf
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    java_bin = args.java_bin
    gatk_jar = args.gatk_jar
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    xmx_mb = args.xmx_mb

    # --- [Output Paths] ---
    if not insert_size_metrics_txt:
    insert_size_metrics_txt = f"{qcResDir}/{SeqID}.insert_size.metrics.txt"
    if not insert_size_hist_pdf:
    insert_size_hist_pdf = f"{qcResDir}/{SeqID}.insert_size.histogram.pdf"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {java_bin} -XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m -jar {gatk_jar} CollectInsertSizeMetrics --INPUT {BamDir}/{SeqID}.{InputSuffix}.bam --OUTPUT {insert_size_metrics_txt} --Histogram_FILE {insert_size_hist_pdf} --REFERENCE_SEQUENCE {ReferenceFasta}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if insert_size_hist_pdf:
        _tgt = os.path.dirname(insert_size_hist_pdf) if os.path.splitext(insert_size_hist_pdf)[1] else insert_size_hist_pdf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if insert_size_metrics_txt:
        _tgt = os.path.dirname(insert_size_metrics_txt) if os.path.splitext(insert_size_metrics_txt)[1] else insert_size_metrics_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()