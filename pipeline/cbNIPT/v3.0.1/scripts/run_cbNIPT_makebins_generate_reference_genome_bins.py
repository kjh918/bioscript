#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = cbNIPT_MakeBins
# VERSION = 1.0.0
# THREADS = 1
# PROFILE = Generate reference genome bins

"""
Tool: cbNIPT_MakeBins (1.0.0)
Profile: Generate reference genome bins
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="cbNIPT_MakeBins Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--ReferenceFasta', required=True, default='', help='참조 유전체 FASTA 파일 경로 (Default: )')
    parser.add_argument('--MappabilityBW', required=True, default='', help='Mappability BigWig 파일 (Default: )')
    parser.add_argument('--OutBinFile', required=True, default='', help='생성된 빈 정보를 저장할 경로 (예: hg38_100kb_bins.bed.gz) (Default: )')
    parser.add_argument('--BinSize', required=False, default='100000', help='Bin 크기 (100kb) (Default: 100000)')
    parser.add_argument('--MinMappability', required=False, default='0.9', help='최소 Mappability 점수 (Default: 0.9)')
    parser.add_argument('--MinGC', required=False, default='0.3', help='최소 GC 함량 필터 (Default: 0.3)')
    parser.add_argument('--MaxGC', required=False, default='0.7', help='최대 GC 함량 필터 (Default: 0.7)')
    parser.add_argument('--IncludeSexChrom', required=False, default='--IncludeSexChrom', help='성염색체 포함 여부 플래그 (제외 시 빈 문자열) (Default: --IncludeSexChrom)')
    parser.add_argument('--python_bin', required=False, default='/storage/home/jhkim/Apps/Python-3.11.13/python', help='Python executable (Default: python)')
    parser.add_argument('--script_path', required=True, default='/storage/home/jhkim/scripts/bioscript/manual/cbnipt_manual_cnv_call/scripts/cli.py', help='Path to the CNV analyzer script (Default: )')
    parser.add_argument('--Threads', required=False, default='1', help='Threads')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    ReferenceFasta = args.ReferenceFasta
    MappabilityBW = args.MappabilityBW
    OutBinFile = args.OutBinFile
    BinSize = args.BinSize
    MinMappability = args.MinMappability
    MinGC = args.MinGC
    MaxGC = args.MaxGC
    IncludeSexChrom = args.IncludeSexChrom
    python_bin = args.python_bin
    script_path = args.script_path

    # --- [Output Paths] ---

    # --- [Command Execution] ---
    cmd = f"{python_bin} {script_path} make-bins --ReferenceFasta {ReferenceFasta} --MappabilityBW {MappabilityBW} --OutBinFile {OutBinFile} --BinSize {BinSize} --MinMappability {MinMappability} --MinGC {MinGC} --MaxGC {MaxGC} {IncludeSexChrom}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if OutBinFile:
        _tgt = os.path.dirname(OutBinFile) if os.path.splitext(OutBinFile)[1] else OutBinFile
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()