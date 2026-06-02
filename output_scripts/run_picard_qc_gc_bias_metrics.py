#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = picard
# VERSION = 3.1.0
# THREADS = 1
# PROFILE = qc_gc_bias_metrics

"""
Tool: picard (3.1.0)
Profile: qc_gc_bias_metrics
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="picard Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for file naming (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM and results (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to the reference genome FASTA file (Default: )')
    parser.add_argument('--gc_bias_metrics_txt', required=False, default='[BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.txt', help='Detailed GC bias metrics for each GC bin (Default: [BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.txt)')
    parser.add_argument('--gc_bias_summary_txt', required=False, default='[BamDir]/[SeqID].[InputSuffix].gc_bias_summary.txt', help='Summary metrics for GC bias (Default: [BamDir]/[SeqID].[InputSuffix].gc_bias_summary.txt)')
    parser.add_argument('--gc_bias_chart_pdf', required=False, default='[BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.pdf', help='PDF chart visualizing normalized coverage by GC content (Default: [BamDir]/[SeqID].[InputSuffix].gc_bias_metrics.pdf)')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of the input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)')
    parser.add_argument('--java_bin', required=False, default='java', help='Path to java executable (Default: java)')
    parser.add_argument('--picard_jar', required=False, default='/storage/apps/bin/picard-3.1.0.jar', help='Path to Picard.jar file (Default: /storage/apps/bin/picard-3.1.0.jar)')
    parser.add_argument('--gc_threads', required=False, default='14', help='Number of threads for Java ParallelGC (Default: 14)')
    parser.add_argument('--xmx_mb', required=False, default='16384', help='Maximum Java heap memory (MB) (Default: 16384)')
    parser.add_argument('--min_gc', required=False, default='0', help='Minimum GC content to include in metrics (Default: 0)')
    parser.add_argument('--max_gc', required=False, default='100', help='Maximum GC content to include in metrics (Default: 100)')
    parser.add_argument('--window_size', required=False, default='100', help='Window size for GC content calculation (Default: 100)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    ReferenceFasta = args.ReferenceFasta
    gc_bias_metrics_txt = args.gc_bias_metrics_txt
    gc_bias_summary_txt = args.gc_bias_summary_txt
    gc_bias_chart_pdf = args.gc_bias_chart_pdf
    InputSuffix = args.InputSuffix
    java_bin = args.java_bin
    picard_jar = args.picard_jar
    gc_threads = args.gc_threads
    xmx_mb = args.xmx_mb
    min_gc = args.min_gc
    max_gc = args.max_gc
    window_size = args.window_size

    # --- [Output Paths] ---
    if not gc_bias_metrics_txt:
    gc_bias_metrics_txt = f"{BamDir}/{SeqID}.{InputSuffix}.gc_bias_metrics.txt"
    if not gc_bias_summary_txt:
    gc_bias_summary_txt = f"{BamDir}/{SeqID}.{InputSuffix}.gc_bias_summary.txt"
    if not gc_bias_chart_pdf:
    gc_bias_chart_pdf = f"{BamDir}/{SeqID}.{InputSuffix}.gc_bias_metrics.pdf"

    # --- [Command Execution] ---
    cmd = f"{java_bin} -XX:ParallelGCThreads={gc_threads} -Xmx{xmx_mb}m -jar {picard_jar} CollectGcBiasMetrics I={BamDir}/{SeqID}.{InputSuffix}.bam O={gc_bias_metrics_txt} S={gc_bias_summary_txt} CHART={gc_bias_chart_pdf} R={ReferenceFasta} MINIMUM_GC={min_gc} MAXIMUM_GC={max_gc} WINDOW_SIZE={window_size}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if gc_bias_chart_pdf:
        _tgt = os.path.dirname(gc_bias_chart_pdf) if os.path.splitext(gc_bias_chart_pdf)[1] else gc_bias_chart_pdf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if gc_bias_summary_txt:
        _tgt = os.path.dirname(gc_bias_summary_txt) if os.path.splitext(gc_bias_summary_txt)[1] else gc_bias_summary_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if gc_bias_metrics_txt:
        _tgt = os.path.dirname(gc_bias_metrics_txt) if os.path.splitext(gc_bias_metrics_txt)[1] else gc_bias_metrics_txt
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()