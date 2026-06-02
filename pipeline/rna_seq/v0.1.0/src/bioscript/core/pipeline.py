import re
import shlex
from pathlib import Path
from typing import Dict, Any, List, Optional
import re
import shlex
import os
from pathlib import Path
from typing import Dict, Any, List

class Task:
    def __init__(self, name: str, runner_path: Path, spec: Dict[str, Any], log_path: Path):
        self.name = name
        self.runner_path = Path(runner_path)
        self.spec = spec
        self.log_path = Path(log_path)
        self.log_path.mkdir(parents=True, exist_ok=True)
        
        # [핵심] 러너 파이썬 파일을 읽고 spec을 주입하여 실제 출력 경로 리스트 생성
        self.outputs = self._parse_outputs_from_runner()
        
        self.exec_type = "/storage/apps/miniconda3/bin/python" if self.runner_path.suffix == ".py" else "bash"

    def _parse_outputs_from_runner(self) -> List[str]:
        """
        러너 스크립트(.py)의 '# --- [Output Paths] ---' 섹션을 분석
        spec 딕셔너리 내용을 바탕으로 {변수}를 실제 값으로 치환함
        """
        if not self.runner_path.exists() or self.runner_path.suffix != ".py":
            return []

        content = self.runner_path.read_text()
        
        # 1. Output 섹션 추출
        pattern = r"# --- \[Output Paths\] ---(.*?)# --- \["
        match = re.search(pattern, content, re.DOTALL)
        if not match:
            return []

        output_section = match.group(1)
        
        # 2. 치환용 변수 맵 준비 (Path 객체는 문자열로 변환)
        # 예: 'qcResDir': PosixPath('/storage/.../fastq_raw') -> '/storage/.../fastq_raw'
        vars_map = {k: str(v) for k, v in self.spec.items()}

        lines = output_section.strip().split('\n')
        resolved_paths = []

        for line in lines:
            if '=' not in line: continue
            
            # 따옴표 내부의 템플릿 추출 (f"{qcResDir}/{SeqID}_R1_fastqc.html")
            path_match = re.search(r'["\'](.*?)["\']', line)
            if path_match:
                template = path_match.group(1)
                
                # 3. f-string 스타일의 {변수}를 spec 값으로 치환
                def replace_var(m):
                    var_name = m.group(1)
                    # spec에 해당 키가 있으면 치환, 없으면 원래 형태 유지
                    return vars_map.get(var_name, f"{{{var_name}}}")
                
                # {qcResDir} -> /storage/home/.../fastq_raw 로 치환됨
                resolved_path = re.sub(r"\{([A-Za-z0-9_]+)\}", replace_var, template)
                resolved_paths.append(resolved_path)

        return resolved_paths

    def get_cmd(self) -> str:
        """spec 딕셔너리를 --key value 형태의 CLI 인자로 변환"""
        args = []
        for k, v in self.spec.items():
            if v is not None:
                # 인자 값 이스케이프 및 문자열화
                args.append(f"--{k} {shlex.quote(str(v))}")
                
                # 'Dir' 파라미터 자동 생성 로직
                if k.endswith('Dir'): 
                    Path(v).mkdir(parents=True, exist_ok=True)
        
        arg_str = " \\\n    ".join(args)
        if args:
            return f"{self.exec_type} {shlex.quote(str(self.runner_path))} \\\n    {arg_str}"
        return f"{self.exec_type} {shlex.quote(str(self.runner_path))}"
    
class Pipeline:
    def __init__(self, raw_dir: str, work_dir: Path, log_dir: Path, executor: Any,max_samples:int,max_threads:int, suffix: str='fastq.gz'):
        self.raw_dir = Path(raw_dir) if raw_dir else None
        self.work_dir = Path(work_dir)
        self.log_dir = Path(log_dir)
        self.executor = executor
        self.max_samples = max_samples
        self.max_threads = max_threads
        self.suffix = suffix
        
        self.sample_tasks: Dict[str, List[Task]] = {}
        self.samples = self._discover_samples(suffix=self.suffix)
        max_samples

    def _discover_samples(self, suffix: str) -> List[str]:
        """RawFastqDir에서 파일 패턴에 따라 샘플 ID 추출"""
        if not self.raw_dir or not self.raw_dir.exists():
            return []
        if suffix == 'bam':
            return sorted(list(set([f.name.split('.')[0] for f in self.raw_dir.glob("*.bam")])))
        return sorted(list(set([f.name.split('_R1')[0] for f in self.raw_dir.glob(f"**/*_R1.{suffix}")])))

    def add_tasks(self, tasks: List[Task]):
        """샘플 ID별로 Task 리스트 그룹화"""
        for t in tasks:
            sid = t.spec.get('SeqID')
            if sid:
                if sid not in self.sample_tasks:
                    self.sample_tasks[sid] = []
                self.sample_tasks[sid].append(t)

    def run(self, run_integrity: bool = True):
        """
        [MODIFIED] 변경 이유: SunGridExecutor 내부의 wait_for_resources()를 호출하여 
        호스트 머신의 활성 슬롯(qstat)을 기준으로 안전하게 자원 제어를 위임함.
        """
        for sid, tasks in self.sample_tasks.items():
            last_jid = None
            all_sample_files = [] 
            master_script_lines = [f"#!/bin/bash", f"# Master Script for {sid}", ""]
            
            # 기완료 여부 확인 (전체 Task가 모두 끝났는지)
            will_run = any(not (task.log_path / ".done").exists() for task in tasks)
            
            if will_run:
                # 1. [핵심] 현재 샘플을 제출하기 전에 Executor에게 자원 확인 요청
                # 샘플 내 최고 스레드 요구량을 기준으로 검사
                sample_peak_threads = max((task.spec.get('Threads', 1) or task.spec.get('threads', 1)) for task in tasks)
                
                print(f"[*] {sid}: Checking SGE resources for {sample_peak_threads} threads...")
                self.executor.wait_for_resources(required_threads=sample_peak_threads)

            # 2. 일반 분석 Task들 순차 제출
            for i, task in enumerate(tasks, start=1):
                all_sample_files.extend(task.outputs)
                job_id = f"{sid}_{i:02d}_{task.name}"
                
                # Skip 로직 (개별 Task)
                if (task.log_path / ".done").exists():
                    master_script_lines.append(f"# Task: {task.name} - SKIPPED")
                    continue
                
                wrapper_p = task.log_path / f"{i:02d}_{task.name}_{sid}.sh"
                self.executor.make_wrapper(task.get_cmd(), wrapper_p, task.log_path)
                
                # SGE 작업 제출 (동일 샘플 내에서는 직전 Task 끝날 때까지 hold_jid 대기)
                last_jid = self.executor.qsub(
                    job_id=job_id,
                    script_path=wrapper_p,
                    log_path=task.log_path,
                    threads=task.spec.get('Threads', 1) or task.spec.get('threads', 1),
                    hold_jid=last_jid 
                )
                
                hold_opt = f"-hold_jid {last_jid}" if last_jid else ""
                master_script_lines.append(f"qsub -N {job_id} {hold_opt} {wrapper_p}")

            # 3. 마지막 Integrity 체크 작업 제출
            if run_integrity and all_sample_files and last_jid:
                self._submit_integrity_task(sid, all_sample_files, last_jid, master_script_lines)

            # 마스터 스크립트 저장 (재현용)
            master_path = self.executor.session_dir / f"run_{sid}.sh"
            master_path.write_text("\n".join(master_script_lines))
            master_path.chmod(0o755)
            
            if will_run:
                print(f"[*] {sid}: Pipeline submitted to SGE. (Last JID: {last_jid})\n")

    def _submit_integrity_task(self, sid, files, hold_jid, master_lines):
        """MD5 체크섬 전용 쉘 스크립트를 생성하고 SGE에 마지막 작업으로 등록"""
        integrity_dir = self.work_dir / sid / "logs" / "99_integrity"
        integrity_dir.mkdir(parents=True, exist_ok=True)
        
        manifest_file = self.work_dir / sid / f"{sid}_final_manifest.md5"
        script_path = integrity_dir / f"run_md5check_{sid}.sh"
        
        # SGE가 실행할 쉘 스크립트 내용 작성
        unique_files = sorted(list(set(files)))
        content = [
            "#!/bin/bash",
            f"# SGE Job for Integrity Check of {sid}",
            f"echo '[$(date)] Starting MD5 check for {sid}' > {integrity_dir}/integrity.log",
            f"rm -f {manifest_file}", # 기존 파일 삭제
            ""
        ]
        
        for f in unique_files:
            # SGE 서버가 실행 시점에 파일이 있는지 체크하고 MD5를 찍음
            content.append(f"if [ -f '{f}' ]; then")
            content.append(f"  md5sum '{f}' >> {manifest_file}")
            content.append(f"else")
            content.append(f"  echo '[MISSING] {f}' >> {integrity_dir}/integrity.log")
            content.append(f"fi")
        
        content.append(f"\necho '[$(date)] All checks finished.' >> {integrity_dir}/integrity.log")
        # 모든 체크가 무사히 끝나면 전체 샘플 완료 표시로 .done 파일 생성 (옵션)
        content.append(f"touch {self.work_dir / sid}/.all_steps_done")

        script_path.write_text("\n".join(content))
        script_path.chmod(0o755)
        
        # SGE에 마지막 작업으로 제출
        final_jid = self.executor.qsub(
            job_id=f"{sid}_integrity",
            script_path=script_path,
            log_path=integrity_dir,
            threads=1,
            hold_jid=hold_jid # 앞선 모든 분석 조브가 끝나야 실행됨
        )
        
        master_lines.append(f"# Final Integrity Task (JID: {final_jid})")
        master_lines.append(f"qsub -N {sid}_integrity -hold_jid {hold_jid} {script_path}\n")