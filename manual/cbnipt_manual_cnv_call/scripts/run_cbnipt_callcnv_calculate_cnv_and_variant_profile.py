#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = cbNIPT_CallCNV
# VERSION = 1.0.0
# THREADS = 4
# PROFILE = Calculate CNV and variant profile

"""
Tool: cbNIPT_CallCNV (1.0.0)
Profile: Calculate CNV and variant profile
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="cbNIPT_CallCNV Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='분석 고유 ID (Default: )')
    parser.add_argument('--BamPath', required=True, default='', help='분석할 단일 BAM 파일 경로 (Default: )')
    parser.add_argument('--VcfFile', required=True, default='', help='')
    parser.add_argument('--AnnotatedBins', required=True, default='', help='make-bins로 생성한 사전 계산 빈 파일(.bed.gz) (Default: )')
    parser.add_argument('--OutDir', required=True, default='', help='결과 저장 디렉토리 (Default: )')
    parser.add_argument('--MinMapQ', required=False, default='20', help='최소 Mapping Quality (Default: 20)')
    parser.add_argument('--LowessFrac', required=False, default='0.2', help='GC 보정 LOWESS 비율 (Default: 0.2)')
    parser.add_argument('--BaselinePloidy', required=False, default='2', help='기준 Ploidy (Default: 2)')
    parser.add_argument('--SmoothWindow', required=False, default='5', help='Rolling median window size (Default: 5)')
    parser.add_argument('--MinDepth', required=False, default='2.0', help='최소 Depth (Default: 2.0)')
    parser.add_argument('--MinCoverage', required=False, default='0.5', help='최소 Coverage 비율 (Default: 0.5)')
    parser.add_argument('--MaskLowerQ', required=False, default='0.01', help='하위 Coverage 필터 (Default: 0.01)')
    parser.add_argument('--MaskUpperQ', required=False, default='0.99', help='상위 Coverage 필터 (Default: 0.99)')
    parser.add_argument('--SegPenalty', required=False, default='10.0', help='Segmentation 민감도 (Default: 10.0)')
    parser.add_argument('--Threads', required=False, default='4', help='사용할 스레드 수 (Default: 4)')
    parser.add_argument('--python_bin', required=False, default='/storage/home/jhkim/Apps/Python-3.11.13/python', help='Python executable (Default: python)')
    parser.add_argument('--script_path', required=True, default='', help='Path to the CNV analyzer script (Default: )')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BamPath = args.BamPath
    VcfFile = args.VcfFile
    AnnotatedBins = args.AnnotatedBins
    OutDir = args.OutDir
    MinMapQ = args.MinMapQ
    LowessFrac = args.LowessFrac
    BaselinePloidy = args.BaselinePloidy
    SmoothWindow = args.SmoothWindow
    MinDepth = args.MinDepth
    MinCoverage = args.MinCoverage
    MaskLowerQ = args.MaskLowerQ
    MaskUpperQ = args.MaskUpperQ
    SegPenalty = args.SegPenalty
    Threads = args.Threads
    python_bin = args.python_bin
    script_path = args.script_path

    # --- [Output Paths] ---

    # --- [Command Execution] ---
    cmd = f"{python_bin} {script_path} call-cnv --SeqID {SeqID} --VcfFile {VcfFile} --BamPath {BamPath} --AnnotatedBins {AnnotatedBins} --OutDir {OutDir} --MinMapQ {MinMapQ} --LowessFrac {LowessFrac} --BaselinePloidy {BaselinePloidy} --SmoothWindow {SmoothWindow} --MinDepth {MinDepth} --MinCoverage {MinCoverage} --MaskLowerQ {MaskLowerQ} --MaskUpperQ {MaskUpperQ} --SegPenalty {SegPenalty} --Threads {Threads}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()