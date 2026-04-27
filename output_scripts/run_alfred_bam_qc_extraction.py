#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = alfred
# VERSION = 0.2.6
# THREADS = 1
# PROFILE = bam_qc_extraction

"""
Tool: alfred (0.2.6)
Profile: bam_qc_extraction
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="alfred Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--TargetBed', required=True, default='', help='Target BED file for WES metrics (Default: )')
    parser.add_argument('--InputSuffix', required=True, default='analysisReady', help='No description (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--alfred_sif', required=False, default='/storage/images/alfred-0.2.6.sif', help='No description (Default: /storage/images/alfred-0.2.6.sif)')
    parser.add_argument('--alfred_bin', required=False, default='/opt/alfred/bin/alfred', help='No description (Default: /opt/alfred/bin/alfred)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    qcResDir = args.qcResDir
    ReferenceFasta = args.ReferenceFasta
    TargetBed = args.TargetBed
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    alfred_sif = args.alfred_sif
    alfred_bin = args.alfred_bin
    bind = args.bind

    # --- [Output Paths] ---
    alfred_raw_tsv = f"{qcResDir}/{SeqID}.alfred.qc.tsv.gz"
    chr_map_stats = f"{qcResDir}/{SeqID}.alfred.chr.map.stats.txt"
    target_coverage = f"{qcResDir}/{SeqID}.alfred.target.coverage.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {alfred_sif} {alfred_bin} qc --reference {ReferenceFasta} --bed {TargetBed} --outfile {alfred_raw_tsv} {BamDir}/{SeqID}.{InputSuffix}.bam && zgrep '^CM' {alfred_raw_tsv} | cut -f 2- > {chr_map_stats} && zgrep '^TC' {alfred_raw_tsv} | cut -f 2- > {target_coverage}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()