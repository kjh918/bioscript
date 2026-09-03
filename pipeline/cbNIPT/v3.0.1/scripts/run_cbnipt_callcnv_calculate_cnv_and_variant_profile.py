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
    parser.add_argument('--AnnotatedBins', required=True, default='', help='make-bins로 생성한 사전 계산 빈 파일(.bed.gz) (Default: )')
    parser.add_argument("--VcfFile", type=str, required=True, help="Path to the merged VCF file (.vcf.gz)")
    parser.add_argument('--OutDir', required=True, default='', help='결과 저장 디렉토리 (Default: )')
    parser.add_argument('--MinMapQ', required=False, default='20', help='최소 Mapping Quality (Default: 20)')
    parser.add_argument('--LowessFrac', required=False, default='0.2', help='GC 보정 LOWESS 비율 (Default: 0.2)')
    parser.add_argument('--BaselinePloidy', required=False, default='2', help='기준 Ploidy (Default: 2)')
    parser.add_argument('--SmoothWindow', required=False, default='5', help='Rolling median window size (Default: 5)')
    parser.add_argument('--MinDepth', required=False, default='2.0', help='최소 Depth (Default: 2.0)')
    parser.add_argument('--MinCoverage', required=False, default='0.2', help='최소 Coverage 비율 (Default: 0.2)')
    parser.add_argument('--MaskLowerQ', required=False, default='0.01', help='하위 Coverage 필터 (Default: 0.01)')
    parser.add_argument('--MaskUpperQ', required=False, default='0.99', help='상위 Coverage 필터 (Default: 0.99)')
    parser.add_argument('--SegPenalty', required=False, default='10.0', help='Segmentation 민감도 (Default: 10.0)')
    parser.add_argument('--Threads', required=False, default='4', help='사용할 스레드 수 (Default: 4)')
    
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
    
    # --- [Output Paths] ---

    # --- [Command Execution] ---
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _package_path = os.path.normpath(
        os.path.join(_script_dir, "..", "src", "cbnipt-cnv-caller")
    )
    _venv_python = os.path.join(_package_path, ".venv", "bin", "python")
    _uv_candidates = [
        os.path.expanduser("~/.local/bin/uv"),
        "/usr/local/bin/uv",
        "/usr/bin/uv",
    ]

    if os.path.exists(_venv_python):
        cmd = (
            f"{_venv_python} -m cbnipt_cnv_caller.cli call-cnv "
            f"--SeqID {SeqID} "
            f"--BamPath {BamPath} "
            f"--VcfFile {VcfFile} "
            f"--AnnotatedBins {AnnotatedBins} "
            f"--OutDir {OutDir} "
            f"--MinMapQ {MinMapQ} "
            f"--LowessFrac {LowessFrac} "
            f"--BaselinePloidy {BaselinePloidy} "
            f"--SmoothWindow {SmoothWindow} "
            f"--MinDepth {MinDepth} "
            f"--MinCoverage {MinCoverage} "
            f"--MaskLowerQ {MaskLowerQ} "
            f"--MaskUpperQ {MaskUpperQ} "
            f"--SegPenalty {SegPenalty} "
            f"--Threads {Threads}"
        )
    else:
        uv_bin = next((p for p in _uv_candidates if os.path.exists(p)), "uv")
        cmd = (
            f"{uv_bin} run --project {_package_path} "
            f"python -m cbnipt_cnv_caller.cli call-cnv "
            f"--SeqID {SeqID} "
            f"--BamPath {BamPath} "
            f"--VcfFile {VcfFile} "
            f"--AnnotatedBins {AnnotatedBins} "
            f"--OutDir {OutDir} "
            f"--MinMapQ {MinMapQ} "
            f"--LowessFrac {LowessFrac} "
            f"--BaselinePloidy {BaselinePloidy} "
            f"--SmoothWindow {SmoothWindow} "
            f"--MinDepth {MinDepth} "
            f"--MinCoverage {MinCoverage} "
            f"--MaskLowerQ {MaskLowerQ} "
            f"--MaskUpperQ {MaskUpperQ} "
            f"--SegPenalty {SegPenalty} "
            f"--Threads {Threads}"
        )

    print(f"[INFO] package_path: {_package_path}")
    print(f"[INFO] venv_python: {_venv_python} ({'found' if os.path.exists(_venv_python) else 'not found, using uv'})")
    if OutDir:
        _tgt = os.path.dirname(OutDir) if os.path.splitext(OutDir)[1] else OutDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()