import shlex
from pathlib import Path
from typing import Dict, Any, List, Optional

class Task:
    def __init__(self, name: str, runner_path: Path, spec: Dict[str, Any], log_path: Path):
        self.name = name
        self.runner_path = Path(runner_path)
        self.spec = spec
        self.log_path = Path(log_path)
        self.log_path.mkdir(parents=True, exist_ok=True)
        # 확장자에 따른 실행 타입 결정
        self.exec_type = "/storage/apps/miniconda3/bin/python" if self.runner_path.suffix == ".py" else "bash"

    def get_cmd(self) -> str:
        """spec 딕셔너리를 --key value 형태의 CLI 인자로 변환 (가독성을 위한 줄바꿈 포함)"""
        args = []
        for k, v in self.spec.items():
            if v is not None:
                # 1. 인자 리스트 구성 (각 인자를 --key 'value' 형태로 한 덩어리씩 저장)
                args.append(f"--{k} {shlex.quote(str(v))}")
                
                # 2. 'Dir'로 끝나는 파라미터는 실행 전 디렉토리 자동 생성
                if k.endswith('Dir'):
                    Path(v).mkdir(parents=True, exist_ok=True)
        
        # 3. 가독성을 위해 백슬래시(\)와 줄바꿈(\n), 들여쓰기(indent)를 넣어 조립
        if args:
            arg_str = " \\\n    ".join(args)
            return f"{self.exec_type} {shlex.quote(str(self.runner_path))} \\\n    {arg_str}"
        else:
            # 인자가 없는 경우
            return f"{self.exec_type} {shlex.quote(str(self.runner_path))}"
        

class Pipeline:
    def __init__(self, raw_dir: str, work_dir: Path, log_dir: Path, executor: Any, suffix: str='fastq'):
        self.raw_dir = Path(raw_dir) if raw_dir else None
        self.work_dir = Path(work_dir)
        self.log_dir = Path(log_dir)
        self.executor = executor
        self.suffix = suffix
        
        # 샘플별 Task 리스트를 저장할 딕셔너리
        self.sample_tasks: Dict[str, List[Task]] = {}
        self.samples = self._discover_samples(suffix=self.suffix)

    def _discover_samples(self, suffix: str) -> List[str]:
        """RawFastqDir에서 R1 파일을 찾아 샘플 ID 추출"""
        if not self.raw_dir or not self.raw_dir.exists():
            return []
        if suffix == 'bam':
            return sorted(list(set([f.name.split('.')[0] for f in self.raw_dir.glob(f"*.bam")])))
        # _R1.fastq.gz 패턴 기준 (필요시 수정)
        return sorted(list(set([f.name.split('_R1')[0] for f in self.raw_dir.glob(f"*_R1.{suffix}")])))
        

    def add_tasks(self, tasks: List[Task]):
        """샘플 ID별로 Task 리스트 그룹화"""
        for t in tasks:
            sid = t.spec.get('SeqID')
            if sid:
                if sid not in self.sample_tasks:
                    self.sample_tasks[sid] = []
                self.sample_tasks[sid].append(t)

    def run(self):
        """샘플별로 각 분석 단계의 logs/에 .sh를 만들고, 중앙에 run_{sid}.sh 생성"""
        for sid, tasks in self.sample_tasks.items():
            last_jid = None
            master_script_lines = [f"#!/bin/bash", f"# Master Script for {sid}", ""]
            
            for i, task in enumerate(tasks, start=1):
                # 1. 작업 디렉토리 (logs의 부모) 및 Job ID
                task_work_dir = task.log_path.parent
                job_id = f"{sid}_{i:02d}_{task.name}"

                # [MODIFIED] .done 파일 존재 시 Skip 로직
                done_flag = task_work_dir / ".done"
                if done_flag.exists():
                    print(f"[-] {job_id}: Already done. Skipping...")
                    master_script_lines.append(f"# Task: {task.name} - SKIPPED (Already Done)")
                    # Skip 시 last_jid는 업데이트하지 않음 (이전 단계 JID 유지 혹은 None)
                    continue
                
                # 2. 개별 Task 스크립트 경로 (요청하신 대로 sid/analysis/logs/ 아래에 생성)
                # 파일명: {order}_{name}_{sid}.sh
                wrapper_p = task.log_path / f"{i:02d}_{task.name}_{sid}.sh"
                
                # 3. 래퍼 스크립트 물리적 생성
                cmd = task.get_cmd()
                self.executor.make_wrapper(cmd, wrapper_p, task_work_dir)
                
                # 4. SGE 제출
                threads = task.spec.get('Threads', 1)
                jid = self.executor.qsub(
                    job_id=job_id,
                    script_path=wrapper_p,
                    log_path=task.log_path,
                    threads=threads,
                    hold_jid=last_jid
                )
                
                # 5. 마스터 로그에 qsub 명령어 기록 (재현용)
                hold_opt = f"-hold_jid {last_jid}" if last_jid else ""
                master_script_lines.append(f"# Task: {task.name} (JID: {jid})")
                master_script_lines.append(
                    f"qsub -N {job_id} -q {self.executor.node} {hold_opt} "
                    f"-pe smp {threads} -V -cwd {wrapper_p}\n"
                )
                
                last_jid = jid

            # 6. 중앙 세션 디렉토리에 샘플 통합 실행 스크립트 저장
            master_path = self.executor.session_dir / f"run_{sid}.sh"
            master_path.write_text("\n".join(master_script_lines))
            master_path.chmod(0o755)
            
            print(f"[*] {sid}: Submitted {len(tasks)} tasks. Master: {master_path}")