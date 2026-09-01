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
            # 지정된 노드(큐)에 던져진 작업들만 가져옵니다.
            out = subprocess.check_output(["qstat", "-q", self.node, "-u", self.user], text=True)
            lines = out.splitlines()[2:] # qstat 헤더 2줄 제외
            
            if not lines:
                return 0, 0

            unique_samples = set()
            total_threads = 0
            
            for line in lines:
                parts = line.split()
                if len(parts) < 6:
                    continue
                
                # 1. 상태(State) 검사: 에러(Eqw)나 삭제 중(d)인 유령 잡 제외
                state = parts[4]
                if 'd' in state or 'E' in state:
                    continue
                
                # 2. 샘플 ID 추출 (Job Name의 첫 부분)
                job_name = parts[2]
                sample_id = job_name.split('_')[0]
                unique_samples.add(sample_id)
                
                # 3. [완벽 파싱] Slots(스레드) 추출 로직
                # SGE에서 slots는 보통 끝에서 1번째 또는 2번째에 있습니다.
                slots = 1
                if not parts[-1].isdigit():
                    # Case 1: 마지막이 '9-12:1' 처럼 Task-ID 묶음인 경우 -> 무조건 그 앞(-2)이 slots
                    if parts[-2].isdigit():
                        slots = int(parts[-2])
                elif parts[-1].isdigit() and parts[-2].isdigit():
                    # Case 2: 마지막도 숫자(Task-ID 1), 그 앞도 숫자(Slots 15)인 경우 -> 뒤에서 2번째가 slots
                    slots = int(parts[-2])
                elif parts[-1].isdigit() and not parts[-2].isdigit():
                    # Case 3: 마지막 숫자(Slots 15) 앞이 큐 이름(all.q@ngsnode1)인 경우 -> 맨 마지막이 slots
                    slots = int(parts[-1])
                    
                total_threads += slots
                
            return len(unique_samples), total_threads
            
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
            f"echo OK > {os.path.dirname(wrapper_path)}/.done"
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