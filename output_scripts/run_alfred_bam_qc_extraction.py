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
    parser.add_argument('--alfred_raw_tsv', required=False, default='[qcResDir]/[SeqID].alfred.qc.tsv.gz', help='No description (Default: [qcResDir]/[SeqID].alfred.qc.tsv.gz)')
    parser.add_argument('--chr_map_stats', required=False, default='[qcResDir]/[SeqID].alfred.chr.map.stats.txt', help='No description (Default: [qcResDir]/[SeqID].alfred.chr.map.stats.txt)')
    parser.add_argument('--target_coverage', required=False, default='[qcResDir]/[SeqID].alfred.target.coverage.txt', help='No description (Default: [qcResDir]/[SeqID].alfred.target.coverage.txt)')
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
    alfred_raw_tsv = args.alfred_raw_tsv
    chr_map_stats = args.chr_map_stats
    target_coverage = args.target_coverage
    singularity_bin = args.singularity_bin
    alfred_sif = args.alfred_sif
    alfred_bin = args.alfred_bin
    bind = args.bind

    # --- [Output Paths] ---
    if not alfred_raw_tsv:
    alfred_raw_tsv = f"{qcResDir}/{SeqID}.alfred.qc.tsv.gz"
    if not chr_map_stats:
    chr_map_stats = f"{qcResDir}/{SeqID}.alfred.chr.map.stats.txt"
    if not target_coverage:
    target_coverage = f"{qcResDir}/{SeqID}.alfred.target.coverage.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {alfred_sif} {alfred_bin} qc --reference {ReferenceFasta} --bed {TargetBed} --outfile {alfred_raw_tsv} {BamDir}/{SeqID}.{InputSuffix}.bam && zgrep '^CM' {alfred_raw_tsv} | cut -f 2- > {chr_map_stats} && zgrep '^TC' {alfred_raw_tsv} | cut -f 2- > {target_coverage}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if chr_map_stats:
        _tgt = os.path.dirname(chr_map_stats) if os.path.splitext(chr_map_stats)[1] else chr_map_stats
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if target_coverage:
        _tgt = os.path.dirname(target_coverage) if os.path.splitext(target_coverage)[1] else target_coverage
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if alfred_raw_tsv:
        _tgt = os.path.dirname(alfred_raw_tsv) if os.path.splitext(alfred_raw_tsv)[1] else alfred_raw_tsv
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if qcResDir:
        _tgt = os.path.dirname(qcResDir) if os.path.splitext(qcResDir)[1] else qcResDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()