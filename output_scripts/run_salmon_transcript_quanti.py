#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = salmon
# VERSION = 1.10.0
# THREADS = 12
# PROFILE = transcript_quanti

"""
Tool: salmon (1.10.0)
Profile: transcript_quanti
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="salmon Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing input BAM files (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='Path to output directory (Default: )')
    parser.add_argument('--SalmonIndex', required=True, default='', help='Path to Salmon index directory (Default: )')
    parser.add_argument('--quant_sf', required=False, default='[OutputDir]/[SeqID]_quant/quant.sf', help='Transcript-level quantification file (Default: [OutputDir]/[SeqID]_quant/quant.sf)')
    parser.add_argument('--quant_log', required=False, default='[OutputDir]/[SeqID]_quant/logs/salmon_quant.log', help='No description (Default: [OutputDir]/[SeqID]_quant/logs/salmon_quant.log)')
    parser.add_argument('--lib_format', required=False, default='[OutputDir]/[SeqID]_quant/lib_format_counts.json', help='Inferred library type info (Default: [OutputDir]/[SeqID]_quant/lib_format_counts.json)')
    parser.add_argument('--Threads', required=False, default='12', help='No description (Default: 12)')
    parser.add_argument('--libType', required=False, default='A', help='Library type (A: Auto-detect) (Default: A)')
    parser.add_argument('--InputSuffix', required=False, default='.Aligned.toTranscriptome.out.bam', help='No description (Default: .Aligned.toTranscriptome.out.bam)')
    parser.add_argument('--extra_args', required=False, default='', help='No description (Default: )')
    parser.add_argument('--samlon_bin', required=False, default='/storage/apps/salmon-1.10.0/bin/salmon', help='No description (Default: /storage/apps/salmon-1.10.0/bin/salmon)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamDir = args.BamDir
    OutputDir = args.OutputDir
    SalmonIndex = args.SalmonIndex
    quant_sf = args.quant_sf
    quant_log = args.quant_log
    lib_format = args.lib_format
    Threads = args.Threads
    libType = args.libType
    InputSuffix = args.InputSuffix
    extra_args = args.extra_args
    samlon_bin = args.samlon_bin

    # --- [Output Paths] ---
    if not quant_sf:
    quant_sf = f"{OutputDir}/{SeqID}_quant/quant.sf"
    if not quant_log:
    quant_log = f"{OutputDir}/{SeqID}_quant/logs/salmon_quant.log"
    if not lib_format:
    lib_format = f"{OutputDir}/{SeqID}_quant/lib_format_counts.json"

    # --- [Command Execution] ---
    cmd = f"{samlon_bin} quant -t {SalmonIndex} --threads {Threads} --libType {libType} -a {BamDir}/{SeqID}{InputSuffix} -o {OutputDir} {extra_args}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if quant_log:
        _tgt = os.path.dirname(quant_log) if os.path.splitext(quant_log)[1] else quant_log
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if OutputDir:
        _tgt = os.path.dirname(OutputDir) if os.path.splitext(OutputDir)[1] else OutputDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if quant_sf:
        _tgt = os.path.dirname(quant_sf) if os.path.splitext(quant_sf)[1] else quant_sf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if lib_format:
        _tgt = os.path.dirname(lib_format) if os.path.splitext(lib_format)[1] else lib_format
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()