#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bbmap_bbsplit
# VERSION = 39.81
# THREADS = 1
# PROFILE = integrated_indexing_and_splitting

"""
Tool: bbmap_bbsplit (39.81)
Profile: integrated_indexing_and_splitting
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bbmap_bbsplit Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sample identifier for naming outputs (Default: )')
    parser.add_argument('--Read1', required=True, default='', help='No description (Default: )')
    parser.add_argument('--Read2', required=True, default='', help='No description (Default: )')
    parser.add_argument('--PrimaryRef', required=True, default='', help='Path to primary reference fasta (Default: )')
    parser.add_argument('--OtherRefsArgs', required=True, default='', help='Formatted string of additional references (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='Directory for results and statistics (Default: )')
    parser.add_argument('--MemoryMB', required=True, default='', help='Java heap memory (e.g., 32768) (Default: )')
    parser.add_argument('--Threads', required=True, default='', help='Number of CPU threads (Default: )')
    parser.add_argument('--BBMapPath', required=True, default='', help='Directory path containing bbsplit.sh (Default: )')
    parser.add_argument('--BaseNamePattern', required=True, default='', help='No description (Default: )')
    parser.add_argument('--ExtraArgs', required=False, default='ambiguous2=all', help='No description (Default: ambiguous2=all)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--bbmap_sif', required=False, default='/storage/home/jhkim/Apps/bbmap_39.81.sif', help='No description (Default: /storage/home/jhkim/Apps/bbmap_39.81.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    Read1 = args.Read1
    Read2 = args.Read2
    PrimaryRef = args.PrimaryRef
    OtherRefsArgs = args.OtherRefsArgs
    OutputDir = args.OutputDir
    MemoryMB = args.MemoryMB
    Threads = args.Threads
    BBMapPath = args.BBMapPath
    BaseNamePattern = args.BaseNamePattern
    ExtraArgs = args.ExtraArgs
    singularity_bin = args.singularity_bin
    bbmap_sif = args.bbmap_sif
    bind = args.bind

    # --- [Output Paths] ---

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {bbmap_sif} bbsplit.sh -Xmx{MemoryMB}M ref_primary={PrimaryRef} {OtherRefsArgs} in1={Read1} in2={Read2} basename={BaseNamePattern} refstats={OutputDir}/{SeqID}.refstats.txt threads={Threads} {ExtraArgs}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(OutputDir) if '.' in os.path.basename(OutputDir) else OutputDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()