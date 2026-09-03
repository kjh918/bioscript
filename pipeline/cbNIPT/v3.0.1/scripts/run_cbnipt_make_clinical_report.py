#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = cbNIPT_KaryotypeReport
# VERSION = 1.0.0
# THREADS = 1
# PROFILE = Generate interactive karyotype HTML report from CNV and NIPT results

"""
Tool: cbNIPT_KaryotypeReport (1.0.0)
Profile: Generate interactive karyotype HTML report from CNV and NIPT results
"""
import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description="cbNIPT_KaryotypeReport Analysis Runner")

    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID',     required=True,  default='', help='분석 고유 ID (Default: )')  # 추가
    parser.add_argument('--cnv',    required=True,  default='', help='CNV 결과 TSV 파일 경로 (Default: )')
    parser.add_argument('--nipt',   required=True,  default='', help='NIPT 결과 TSV 파일 경로 (Default: )')
    parser.add_argument('--report', required=True,  default='', help='출력 리포트 경로 (확장자 제외, e.g. output/SAMPLE_001) (Default: )')
    parser.add_argument('--chrom',  required=False, default='1',    help='초기 표시 염색체 번호 (Default: 1)')
    parser.add_argument('--host',   required=False, default='0.0.0.0', help='서버 호스트 (Default: 0.0.0.0)')
    parser.add_argument('--port',   required=False, default='8050', help='서버 포트 (Default: 8050)')
    parser.add_argument('--no-server', action='store_true',
                        help='리포트 파일만 생성하고 서버를 시작하지 않음 (Default: False)')
    parser.add_argument('--Threads', required=False, default='1', help='스레드 수 (Default: 1)')

    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    CnvFile    = args.cnv
    nipt   = args.nipt
    ReportPath = args.report
    Chrom      = args.chrom
    Host       = args.host
    Port       = args.port
    NoServer   = args.no_server
    Threads  = args.Threads  # 추가

    # --- [Output Paths] ---
    _report_dir = os.path.dirname(ReportPath) if os.path.dirname(ReportPath) else "."

    # --- [Command Execution] ---
    _script_dir   = os.path.dirname(os.path.abspath(__file__))
    _package_path = os.path.normpath(
        os.path.join(_script_dir, "..", "src", "cbnipt-report")
    )
    _venv_python  = os.path.join(_package_path, ".venv", "bin", "python")
    _uv_candidates = [
        os.path.expanduser("~/.local/bin/uv"),
        "/usr/local/bin/uv",
        "/usr/bin/uv",
    ]

    # 공통 플래그 조립
    _common_flags = (
        f"--cnv    {CnvFile} "
        f"--nipt   {nipt} "
        f"--report {ReportPath} "
        f"--chrom  {Chrom} "
        f"--host   {Host} "
        f"--port   {Port} "
        + ("--no-server " if NoServer else "")
    )

    if os.path.exists(_venv_python):
        cmd = (
            f"{_venv_python} -m karyotype_viewer "
            + _common_flags
        )
    else:
        uv_bin = next((p for p in _uv_candidates if os.path.exists(p)), "uv")
        cmd = (
            f"{uv_bin} run --project {_package_path} "
            f"python -m karyotype_viewer "
            + _common_flags
        )

    print(f"[INFO] package_path : {_package_path}")
    print(f"[INFO] venv_python  : {_venv_python} "
          f"({'found' if os.path.exists(_venv_python) else 'not found, using uv'})")
    print(f"[INFO] report_dir   : {_report_dir}")

    if _report_dir and _report_dir != ".":
        os.makedirs(_report_dir, exist_ok=True)

    subprocess.run(cmd, shell=True, check=True)


if __name__ == "__main__":
    main()