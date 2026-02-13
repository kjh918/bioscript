from __future__ import annotations
import argparse
import sys
from pathlib import Path
from bioscript.generator.factory import BioPipelineFactory

def main() -> None:
    # 1. 메인 파서 설정
    ap = argparse.ArgumentParser(
        prog="bioscript",
        description="BioScript: Standalone Script Generator for Bioinformatics Pipelines"
    )
    sub = ap.add_subparsers(dest="cmd", required=True, help="Available commands")

    # 2. 'make' 커맨드 설정 (스크립트 생성)
    p_make = sub.add_parser("make", help="Generate standalone scripts (Py, Sh, Nf) from YAML config")
    p_make.add_argument("-c", "--config", required=True, help="Path to the tool configuration YAML file")
    p_make.add_argument("-o", "--outdir", help="Override the output directory (default: ./output_scripts)")

    # 3. 'run' 커맨드 (추후 배포판 실행기를 위해 미리 예약 가능)
    # p_run = sub.add_parser("run", help="Execute generated scripts using DistributedPipeline")

    args = ap.parse_args()

    # --- 실행 로직 ---
    if args.cmd == "make":
        config_path = Path(args.config).resolve()
        
        if not config_path.exists():
            print(f"Error: Config file not found at {config_path}")
            sys.exit(1)

        # 팩토리 초기화
        factory = BioPipelineFactory(config_path)
        
        # 출력 디렉토리 결정 (CLI 옵션 우선 -> 없으면 기본값)
        # YAML 내부의 output_dir 설정을 가져오고 싶다면 factory 내부 로직 활용 가능
        output_dir = Path(args.outdir).resolve() if args.outdir else Path("./output_scripts")

        print(f"[*] Starting build process for: {config_path.name}")
        try:
            factory.build_all(output_dir)
            print(f"\n[✔] Successfully generated scripts in: {output_dir}")
        except Exception as e:
            print(f"Error during build: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()