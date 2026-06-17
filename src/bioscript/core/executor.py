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

    def _get_node_status(self) -> tuple[int, int]:
        """해당 노드(Queue)의 (실제 사용 중인 슬롯, 전체 슬롯) 반환"""
        used, tot = 0, self.max_threads
        try:
            # -f 옵션으로 노드 자체의 리소스 현황을 강제 출력
            out = subprocess.check_output(["qstat", "-f", "-q", self.node], text=True)
            for line in out.splitlines()[2:]:
                if self.node in line:
                    parts = line.split()
                    # 포맷: all.q@ngsnode1  BIP   0/20/40  ...
                    for p in parts:
                        if p.count('/') == 2:
                            resv, u, t = p.split('/')
                            used = int(u)
                            # 유저 설정값과 노드 물리 한계값 중 더 보수적인 것을 총량으로 잡음
                            tot = min(int(t), self.max_threads)
                            break
                    break
        except Exception:
            pass
        return used, tot

    def _get_my_usage(self) -> tuple[int, int]:
        """이 스크립트가 쏜 작업 중 (활성 샘플 수, qw 대기열에 쌓은 스레드 수) 반환"""
        my_samples = set()
        my_qw_threads = 0
        
        if not self.active_jids:
            return 0, 0
            
        try:
            out = subprocess.check_output(["qstat", "-u", self.user], text=True)
            lines = out.splitlines()[2:]
            
            for line in lines:
                parts = line.split()
                if len(parts) < 6: 
                    continue
                
                jid = parts[0]
                # [완벽 차단] 이 스크립트에서 qsub한 내역이 없으면 무조건 스킵
                if jid not in self.active_jids:
                    continue
                    
                state = parts[4]
                if 'd' in state or 'E' in state:
                    continue
                
                job_name = parts[2]
                sample_id = job_name.split('_')[0]
                my_samples.add(sample_id)
                
                # SGE에서 대기 중(qw)인 작업은 아직 노드 사용량(used)에 잡히지 않으므로
                # 우리가 따로 합산해서 예상 부하를 계산해야 합니다.
                if 'qw' in state:
                    slots = 1
                    if not parts[-1].isdigit() and parts[-2].isdigit():
                        slots = int(parts[-2])
                    elif parts[-1].isdigit() and parts[-2].isdigit():
                        slots = int(parts[-2])
                    elif parts[-1].isdigit() and not parts[-2].isdigit():
                        slots = int(parts[-1])
                    my_qw_threads += slots
                    
        except Exception:
            pass
            
        return len(my_samples), my_qw_threads

    def wait_for_resources(self, required_threads: int):
        """가용 자원이 생길 때까지 Polling하며 대기"""
        while True:
            # 1. 노드의 실제 사용량 (다른 사람 작업 포함)
            node_used, node_tot = self._get_node_status()
            # 2. 내 파이프라인의 진행 상태 (내 샘플 수, 내가 대기시킨 스레드)
            my_samples, my_qw_threads = self._get_my_usage()
            
            # 샘플 수 제한 체크
            sample_ok = my_samples < self.max_samples
            
            # 노드 빈 공간 체크 (현재 사용량 + 내 대기열 + 지금 넣을 작업) <= 총량
            projected_used = node_used + my_qw_threads + required_threads
            threads_ok = projected_used <= node_tot
            
            if sample_ok and threads_ok:
                break
                
            print(f"[*] [Waiting on {self.node}] Samples: {my_samples}/{self.max_samples} | "
                  f"Node Load: {node_used}(used) + {my_qw_threads}(my qw) + {required_threads}(new) > {node_tot}. Sleep 30s...")
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