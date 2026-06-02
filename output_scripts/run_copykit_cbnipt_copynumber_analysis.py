#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = CopyKit
# VERSION = 1.0.0
# THREADS = 1
# PROFILE = cbNIPT_CopyNumber_Analysis

"""
Tool: CopyKit (1.0.0)
Profile: cbNIPT_CopyNumber_Analysis
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="CopyKit Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Original sequence identifier from mapping stage (Default: )')
    parser.add_argument('--NGS_DataBaseDir', required=True, default='', help='Directory containing the processed BAM files (Default: )')
    parser.add_argument('--ResultBaseDir', required=True, default='', help='Base directory where all sample results are stored (Default: )')
    parser.add_argument('--AnalysisRunDir', required=False, default='[ResultBaseDir]/[SeqID]', help='Target directory for this specific sample analysis (Default: [ResultBaseDir]/[SeqID])')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)')
    parser.add_argument('--mkdir_bin', required=False, default='mkdir', help='Path to mkdir executable (Default: mkdir)')
    parser.add_argument('--rscript_bin', required=False, default='Rscript', help='Path to Rscript executable (Default: Rscript)')
    parser.add_argument('--Rscript_path', required=False, default='/storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R', help='Path to the CopyKit analysis R script (Default: /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_Run.CopyKit.Analysis.R)')
    parser.add_argument('--Threads', required=False, default='1', help='Number of threads for parallel processing (Default: 1)')
    parser.add_argument('--BinSize', required=False, default='220kb', help='Bin size for copy number estimation (e.g., 220kb, 500kb) (Default: 220kb)')
    parser.add_argument('--GenomeVersion', required=False, default='hg38', help='Reference genome version (hg38/hg19) (Default: hg38)')
    parser.add_argument('--SamplePloidy', required=False, default='2', help='Expected sample ploidy (Default: 2 for human) (Default: 2)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    NGS_DataBaseDir = args.NGS_DataBaseDir
    ResultBaseDir = args.ResultBaseDir
    AnalysisRunDir = args.AnalysisRunDir
    InputSuffix = args.InputSuffix
    mkdir_bin = args.mkdir_bin
    rscript_bin = args.rscript_bin
    Rscript_path = args.Rscript_path
    Threads = args.Threads
    BinSize = args.BinSize
    GenomeVersion = args.GenomeVersion
    SamplePloidy = args.SamplePloidy

    # --- [Output Paths] ---
    if not AnalysisRunDir:
    AnalysisRunDir = f"{ResultBaseDir}/{SeqID}"

    # --- [Command Execution] ---
    cmd = f"{mkdir_bin} -p {AnalysisRunDir} && ln -Tsf {NGS_DataBaseDir}/{SeqID}.{InputSuffix}.bam {AnalysisRunDir}/{SeqID}.bam && ln -Tsf {NGS_DataBaseDir}/{SeqID}.{InputSuffix}.bam.bai {AnalysisRunDir}/{SeqID}.bam.bai && {rscript_bin} {Rscript_path}  --SeqID {SeqID}  --AnalysisRunDir {AnalysisRunDir}  --BinSize {BinSize}  --Ploidy {SamplePloidy}  --GenomeVersion {GenomeVersion}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if NGS_DataBaseDir:
        _tgt = os.path.dirname(NGS_DataBaseDir) if os.path.splitext(NGS_DataBaseDir)[1] else NGS_DataBaseDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if ResultBaseDir:
        _tgt = os.path.dirname(ResultBaseDir) if os.path.splitext(ResultBaseDir)[1] else ResultBaseDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if AnalysisRunDir:
        _tgt = os.path.dirname(AnalysisRunDir) if os.path.splitext(AnalysisRunDir)[1] else AnalysisRunDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()