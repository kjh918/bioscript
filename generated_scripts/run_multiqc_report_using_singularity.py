#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = multiqc
# VERSION = 1.16
# THREADS = 1
# PROFILE = report_using_singularity

"""
Tool: multiqc (1.16)
Profile: report_using_singularity
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="multiqc Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier used for the report filename (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='Directory containing all QC data and target output directory (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--multiqc_bin', required=False, default='multiqc', help='MultiQC binary name inside container (Default: multiqc)')
    parser.add_argument('--sif', required=False, default='/storage/images/multiqc-1.16.sif', help='MultiQC SIF image path (Default: /storage/images/multiqc-1.16.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--mqc_config', required=False, default='/storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml', help='Path to custom MultiQC config file (Default: /storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml)')
    parser.add_argument('--mqc_filename', required=False, default='[SeqID].QC.Results', help='Base name for the output report file (Default: [SeqID].QC.Results)')
    parser.add_argument('--mqc_args', required=False, default='--force --data-dir', help='Additional MultiQC flags (Default: --force --data-dir)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    qcResDir = args.qcResDir
    singularity_bin = args.singularity_bin
    multiqc_bin = args.multiqc_bin
    sif = args.sif
    bind = args.bind
    mqc_config = args.mqc_config
    mqc_filename = args.mqc_filename
    mqc_args = args.mqc_args

    # --- [Output Paths] ---
    report_html = f"{qcResDir}/{SeqID}.QC.Results.html"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {multiqc_bin} {mqc_args} --filename {mqc_filename} --outdir {qcResDir} --config {mqc_config} {qcResDir}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()