#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bedtools_genomecov_to_bigwig
# VERSION = v2.27.1
# THREADS = 4
# PROFILE = genome_coverage_and_signal_generation

"""
Tool: bedtools_genomecov_to_bigwig (v2.27.1)
Profile: genome_coverage_and_signal_generation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bedtools_genomecov_to_bigwig Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Unique sample identifier used to generate the output file names. (Default: )')
    parser.add_argument('--InputFile', required=True, default='', help='Path to the input alignment (BAM) or interval file (BED/GFF/VCF) to compute coverage from. (Default: )')
    parser.add_argument('--GenomeSizes', required=True, default='', help='Path to the tab-delimited chromosome sizes reference file. (Default: )')
    parser.add_argument('--Extension', required=True, default='', help='Intermediate file extension reflecting bedGraph data structure (must be \'bedgraph\' or \'bg\'). (Default: )')
    parser.add_argument('--OutputDir', required=True, default='', help='Target directory where the generated genome coverage (bedGraph and BigWig) files will be saved. (Default: )')
    parser.add_argument('--Threads', required=False, default='4', help='No description (Default: 4)')
    parser.add_argument('--InputFlag', required=False, default='-ibam', help='The conditional input flag choice: Use \'-ibam\' for BAM files or \'-i\' for interval files. (Default: -ibam)')
    parser.add_argument('--GenomeSizesFlag', required=False, default='', help='The chromosome layout configuration string for bedtools: Pass \'-g [GenomeSizes]\' for interval files or \'\' for BAM files. (Default: )')
    parser.add_argument('--ScaleArgs', required=False, default='-bg', help='Pre-calculated scale arguments (e.g., \'-scale 0.05 -bg\'). Note: \'-bg\' is required to match bedGraph format specifications. (Default: -bg)')
    parser.add_argument('--SortCmd', required=False, default='| LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n', help='Downstream shell sort pipeline instruction. Pass \'\' if sorting is disabled. (Default: | LC_ALL=C sort --buffer-size=8G -k1,1 -k2,2n)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Absolute path to the singularity/apptainer runtime executable binary. (Default: singularity)')
    parser.add_argument('--bedtools_sif', required=False, default='/storage/images/bedtools-2.27.1.sif', help='Absolute path to the Bedtools Singularity Image File (.sif). (Default: /storage/images/bedtools-2.27.1.sif)')
    parser.add_argument('--bigwig_sif', required=False, default='/storage/images/ucsc-bedgraphtobigwig-445.sif', help='Absolute path to the Singularity Image File containing the bedGraphToBigWig UCSC utility. (Default: /storage/images/ucsc-bedgraphtobigwig-445.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Comma-separated host-to-container file system mount paths. (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    InputFile = args.InputFile
    GenomeSizes = args.GenomeSizes
    Extension = args.Extension
    OutputDir = args.OutputDir
    Threads = args.Threads
    InputFlag = args.InputFlag
    GenomeSizesFlag = args.GenomeSizesFlag
    ScaleArgs = args.ScaleArgs
    SortCmd = args.SortCmd
    singularity_bin = args.singularity_bin
    bedtools_sif = args.bedtools_sif
    bigwig_sif = args.bigwig_sif
    bind = args.bind

    # --- [Output Paths] ---

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {bedtools_sif} bedtools genomecov {InputFlag} {InputFile} {GenomeSizesFlag} {ScaleArgs} {SortCmd} --parallel={Threads} > {OutputDir}/{SeqID}.{Extension} && {singularity_bin} exec -B {bind} {bigwig_sif} bedGraphToBigWig {OutputDir}/{SeqID}.{Extension} {GenomeSizes} {OutputDir}/{SeqID}.bw"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if OutputDir:
        _tgt = os.path.dirname(OutputDir) if os.path.splitext(OutputDir)[1] else OutputDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()