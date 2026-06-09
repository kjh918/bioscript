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
    def __init__(self, raw_dir: str, work_dir: Path, log_dir: Path, executor: Any, suffix: str='fastq.gz'):
        self.raw_dir = Path(raw_dir) if raw_dir else None
        self.work_dir = Path(work_dir)
        self.log_dir = Path(log_dir)
        self.executor = executor
        self.suffix = suffix
        
        # [MODIFIED] Global Task(샘플 무관)와 Sample Task를 분리하여 저장
        self.global_tasks: List[Task] = []
        self.sample_tasks: Dict[str, List[Task]] = {}
        self.samples = self._discover_samples(suffix=self.suffix)

    def _discover_samples(self, suffix: str) -> List[str]:
        # (기존 코드 동일)
        if not self.raw_dir or not self.raw_dir.exists():
            return []
        if suffix == 'bam':
            return sorted(list(set([f.name.split('.')[0] for f in self.raw_dir.glob("*.bam")])))
        return sorted(list(set([f.name.split('_R1')[0] for f in self.raw_dir.glob(f"**/*_R1.{suffix}")])))

    def add_tasks(self, tasks: List[Task]):
        """샘플 ID별로 Task 리스트 그룹화 및 Global Task 분리"""
        for t in tasks:
            sid = t.spec.get('SeqID')
            if sid:
                # SeqID가 있으면 샘플 종속 Task
                if sid not in self.sample_tasks:
                    self.sample_tasks[sid] = []
                self.sample_tasks[sid].append(t)
            else:
                # [MODIFIED] SeqID가 없으면 1회성 Global Task로 분류
                self.global_tasks.append(t)

    def run(self, run_integrity: bool = True):
        global_last_jid = None
        
        # =========================================================
        # 1. Global Tasks (Reference Indexing, DB Setup 등) 먼저 처리
        # =========================================================
        if self.global_tasks:
            print("[*] Submitting Global/Reference Tasks...")
            master_global_lines = ["#!/bin/bash", "# Master Script for Global Tasks", ""]
            
            for i, task in enumerate(self.global_tasks, start=1):
                job_id = f"GLOBAL_{i:02d}_{task.name}"
                
                if (task.log_path / ".done").exists():
                    master_global_lines.append(f"# Global Task: {task.name} - SKIPPED")
                    continue
                
                wrapper_p = task.log_path / f"{i:02d}_{task.name}_global.sh"
                self.executor.make_wrapper(task.get_cmd(), wrapper_p, task.log_path)
                
                # Global Task들을 순차적으로(chaining) 제출
                global_last_jid = self.executor.qsub(
                    job_id=job_id,
                    script_path=wrapper_p,
                    log_path=task.log_path,
                    threads=task.spec.get('Threads') or task.spec.get('threads') or 1,
                    hold_jid=global_last_jid
                )
                
                hold_opt = f"-hold_jid {global_last_jid}" if global_last_jid else ""
                master_global_lines.append(f"qsub -N {job_id} {hold_opt} {wrapper_p}")
            
            # Global Task용 마스터 스크립트 별도 저장
            global_master_path = self.executor.session_dir / "run_00_GLOBAL_TASKS.sh"
            global_master_path.write_text("\n".join(master_global_lines))
            global_master_path.chmod(0o755)
            print(f"[*] Global tasks submitted. (Last JID: {global_last_jid})")

        # =========================================================
        # 2. Sample Tasks 처리 (Global Task 완료 시점부터 시작)
        # =========================================================
        for sid, tasks in self.sample_tasks.items():
            print(f"[*] Processing Sample: {sid}")
            
            # [MODIFIED] 핵심: 샘플의 첫 작업은 Global Task가 끝날 때까지 대기(-hold_jid)해야 함
            last_jid = global_last_jid 
            
            all_sample_files = [] 
            master_script_lines = [f"#!/bin/bash", f"# Master Script for {sid}", ""]
            
            for i, task in enumerate(tasks, start=1):
                all_sample_files.extend(task.outputs)
                job_id = f"{sid}_{i:02d}_{task.name}"
                
                if (task.log_path / ".done").exists():
                    master_script_lines.append(f"# Task: {task.name} - SKIPPED")
                    continue
                
                wrapper_p = task.log_path / f"{i:02d}_{task.name}_{sid}.sh"
                self.executor.make_wrapper(task.get_cmd(), wrapper_p, task.log_path)
                
                # last_jid가 초기엔 global_last_jid를 가리키므로 첫 샘플 작업이 안전하게 대기함
                last_jid = self.executor.qsub(
                    job_id=job_id,
                    script_path=wrapper_p,
                    log_path=task.log_path,
                    threads=task.spec.get('Threads') or task.spec.get('threads') or 1,
                    hold_jid=last_jid 
                )
                
                hold_opt = f"-hold_jid {last_jid}" if last_jid else ""
                master_script_lines.append(f"qsub -N {job_id} {hold_opt} {wrapper_p}")

            # 무결성 검사 제출 (기존과 동일)
            if run_integrity and all_sample_files and last_jid:
                self._submit_integrity_task(sid, all_sample_files, last_jid, master_script_lines)

            # 샘플별 마스터 스크립트 저장
            master_path = self.executor.session_dir / f"run_{sid}.sh"
            master_path.write_text("\n".join(master_script_lines))
            master_path.chmod(0o755)
            
            print(f"[*] {sid}: Pipeline submitted to SGE. (Last JID: {last_jid})")

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