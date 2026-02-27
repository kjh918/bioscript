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
    p_pipe = sub.add_parser("pipe", help="Generate standalone scripts (Py, Sh, Nf) from YAML config")
    p_pipe.add_argument("-c", "--config", required=True, help="Path to the tool configuration YAML file")
    p_pipe.add_argument("-o", "--outdir", help="Override the output directory (default: ./output_scripts)")

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
# bioscript/cli/main.py (일부 발췌)

    elif args.cmd == "pipe":
        config_path = Path(args.config).resolve()
        if not config_path.exists():
            print(f"Error: Pipeline config file not found at {config_path}")
            sys.exit(1)

        print(f"[*] Initializing Pipeline Orchestrator with: {config_path.name}")
        
        # 1. Pipeline Config 로드 (전역 리소스 및 경로 설정)
        import yaml
        with open(config_path, 'r') as f:
            pipe_cfg = yaml.safe_load(f)

        # 2. 필수 엔진 임포트
        from bioscript.core.executor import SunGridExecutor
        from bioscript.core.pipeline import Pipeline, Task

        # 3. Pipeline 엔진 초기화
        work_dir = Path(pipe_cfg['project']['work_dir'])
        scripts_dir = Path(pipe_cfg['project']['scripts_dir'])
        
        pipe = Pipeline(
            raw_dir = Path(pipe_cfg['project']['raw_dir']),
            work_dir = work_dir,
            log_dir = work_dir / "logs",
            executor = SunGridExecutor(
                user = pipe_cfg['executor']['user'],
                node = pipe_cfg['executor']['node'],
                log_root = work_dir / "logs",
                max_samples = pipe_cfg['executor'].get('max_samples', 5),
                max_threads = pipe_cfg['executor'].get('max_threads', 40)
            )
        )

        # 4. 샘플별 Task 할당 루프
        # YAML에 정의된 순서대로 Task 리스트를 생성하여 연결합니다.
        for sid in pipe.samples:
            sample_tasks = []
            
            # YAML의 'workflow' 섹션에 정의된 도구 순서대로 순회
            for step in pipe_cfg['workflow']:
                tool_name = step['tool']
                # 생성된 스크립트 경로 (예: run_fastp.py)
                runner = scripts_dir / f"run_{tool_name}.py"
                
                # 각 단계별 특화 인자(spec) 구성
                spec = { 'SeqID': sid }
                spec.update(step.get('params', {}))
                
                # 경로 기반 인자들은 샘플 ID(sid)를 포함하여 동적 생성
                # 예: work_dir/sid/bam
                for k, v in step.get('io', {}).items():
                    spec[k] = str(Path(v.replace("[sid]", sid)).absolute())

                sample_tasks.append(
                    Task(
                        name = tool_name,
                        runner_path = runner,
                        log_path = work_dir / sid / "logs" / f"step_{tool_name}",
                        spec = spec
                    )
                )
            
            pipe.add_tasks(sample_tasks)

        # 5. 실행 (SGE hold_jid 자동 관리)
        print(f"[*] Submitting {len(pipe.samples)} samples to SGE Cluster...")
        pipe.run()

    else:
        pass 
if __name__ == "__main__":
    main()