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
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the sorted input BAM files (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory where duplication metrics will be saved (Default: )')
    parser.add_argument('--java_bin', required=False, default='java', help='Path to java executable (Default: java)')
    parser.add_argument('--picard_jar', required=False, default='/storage/apps/bin/picard.jar', help='Path to picard.jar file (Default: /storage/apps/bin/picard.jar)')
    parser.add_argument('--Threads', required=False, default='14', help='Number of parallel GC threads (-XX:ParallelGCThreads) (Default: 14)')
    parser.add_argument('--Memory', required=False, default='16384m', help='Maximum Java heap size (-Xmx) (Default: 16384m)')
    parser.add_argument('--TmpDir', required=False, default='/tmp', help='Directory for Picard temporary files (Default: /tmp)')
    parser.add_argument('--create_index', required=False, default='true', help='Automatically create BAM index (.bai) (Default: true)')
    parser.add_argument('--remove_duplicates', required=False, default='true', help='If true, removes duplicates completely rather than just marking them (Default: true)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    java_bin = args.java_bin
    picard_jar = args.picard_jar
    Threads = args.Threads
    Memory = args.Memory
    TmpDir = args.TmpDir
    create_index = args.create_index
    remove_duplicates = args.remove_duplicates

    # --- [Output Paths] ---
    out_bam = f"{BamDir}/{SeqID}.sorted.dedup.bam"
    metrics = f"{qcResDir}/{SeqID}.mark.duplicates.metrics.txt"

    # --- [Command Execution] ---
    cmd = f"{java_bin} -XX:ParallelGCThreads={Threads} -Xmx{Memory} -jar {picard_jar} MarkDuplicates  INPUT={BamDir}/{SeqID}.sorted.bam  OUTPUT={out_bam}  METRICS_FILE={metrics}  CREATE_INDEX={create_index}  REMOVE_DUPLICATES={remove_duplicates}  TMP_DIR={TmpDir}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    os.makedirs(os.path.dirname(TmpDir) if '.' in os.path.basename(TmpDir) else TmpDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()