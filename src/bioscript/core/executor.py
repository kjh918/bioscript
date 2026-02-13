import os
import subprocess
import shlex
import time
from pathlib import Path

class SunGridExecutor:
    def __init__(self, user: str, node: str, log_root: Path, max_samples: int = 5, max_threads: int = 40):
        self.user = user
        self.node = node
        self.log_root = Path(log_root)
        self.max_samples = max_samples
        self.max_threads = max_threads
        self.run_id = time.strftime("%Y%m%d_%H%M%S")
        self.session_dir = self.log_root / self.run_id
        self.session_dir.mkdir(parents=True, exist_ok=True)

    def _get_active_usage(self) -> tuple[int, int]:
        """현재 SGE 상의 (진행 중인 샘플 수, 사용 중인 총 스레드 수) 반환"""
        try:
            # qstat 출력을 가져와서 파싱
            out = subprocess.check_output(["qstat", "-u", self.user], text=True)
            lines = out.splitlines()[2:] # 헤더 제외
            
            active_jids = [line.split()[0] for line in lines]
            if not active_jids:
                return 0, 0

            # 상세 정보(스레드 수) 확인을 위해 qstat -j 실행
            total_threads = 0
            # Job ID 중복 제거 (동일 샘플 내 여러 작업이 보일 수 있음)
            unique_sample_count = len(set([line.split()[2].split('_')[0] for line in lines]))
            
            # 실무적으로는 qstat 출력의 'slots' 컬럼을 합산
            for line in lines:
                parts = line.split()
                # SGE qstat 출력 컬럼: job-ID, prior, name, user, state, submit/start at, queue, slots, ja-task-ID
                slots = int(parts[-2]) if parts[-2].isdigit() else 1
                total_threads += slots
                
            return unique_sample_count, total_threads
        except (subprocess.CalledProcessError, IndexError):
            return 0, 0

    def wait_for_resources(self, required_threads: int):
        """가용 자원이 생길 때까지 Polling하며 대기"""
        while True:
            current_samples, current_threads = self._get_active_usage()
            
            # 1. 샘플 수 제한 체크
            sample_ok = current_samples < self.max_samples
            # 2. 총 스레드 수 제한 체크 (현재 스레드 + 추가될 스레드)
            threads_ok = (current_threads + required_threads) <= self.max_threads
            
            if sample_ok and threads_ok:
                break
                
            print(f"[*] [Waiting] Resources Full: Samples({current_samples}/{self.max_samples}), Threads({current_threads}/{self.max_threads}). Sleep 30s...")
            time.sleep(30)

    def make_wrapper(self, cmd: str, wrapper_path: Path, work_dir: Path) -> Path:
        wrapper_path.parent.mkdir(parents=True, exist_ok=True)
        script = [
            "#!/bin/bash",
            "#$ -S /bin/bash",
            "set -e",
            cmd,
            "echo OK > .done"
        ]
        wrapper_path.write_text("\n".join(script))
        wrapper_path.chmod(0o755)
        return wrapper_path

    def qsub(self, job_id: str, script_path: Path, log_path: Path, threads: int, hold_jid: str = None) -> str:
        log_path.mkdir(parents=True, exist_ok=True)

        args = [
            "qsub", "-N", job_id, "-q", self.node,
            "-pe", "smp", str(threads),
            "-o", str(log_path).replace("sh","stdout"), "-e", str(log_path).replace("sh","stderr"),
            "-V", "-cwd"
        ]
        if hold_jid:
            args += ["-hold_jid", hold_jid]
        
        args.append(str(script_path))
        out = subprocess.check_output(args, text=True)
        return out.split()[2]