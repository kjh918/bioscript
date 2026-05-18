import os
import time
import subprocess
from pathlib import Path
from typing import Any, Tuple

class SunGridExecutor:
    # [MODIFIED] Removed all default parameter values to enforce Rule 0. Every boundary must be explicitly injected.
    def __init__(self, user: str, node: str, log_root: Path, max_samples: int, max_threads: int):
        # [MODIFIED] Added strict type and value validation blocks to catch misconfigurations before cluster engagement.
        if not user or not isinstance(user, str):
            raise ValueError("[STRICT_ERROR] Parameter 'user' must be a non-empty string.")
        if not node or not isinstance(node, str):
            raise ValueError("[STRICT_ERROR] Parameter 'node' must be a non-empty string.")
        if not log_root or not isinstance(log_root, Path):
            raise ValueError("[STRICT_ERROR] Parameter 'log_root' must be a valid Path object.")
        if max_samples is None or not isinstance(max_samples, int) or max_samples <= 0:
            raise ValueError("[STRICT_ERROR] Parameter 'max_samples' must be an integer greater than 0.")
        if max_threads is None or not isinstance(max_threads, int) or max_threads <= 0:
            raise ValueError("[STRICT_ERROR] Parameter 'max_threads' must be an integer greater than 0.")

        self.user = user
        self.node = node
        self.log_root = log_root
        self.max_samples = max_samples
        self.max_threads = max_threads
        self.run_id = time.strftime("%Y%m%d_%H%M%S")
        self.session_dir = self.log_root / self.run_id
        self.session_dir.mkdir(parents=True, exist_ok=True)

    def _get_active_usage(self) -> Tuple[int, int]:
        """Calculates current SGE resource utilization: (active_sample_count, active_slots_count)"""
        try:
            out = subprocess.check_output(["qstat", "-u", self.user], text=True)
            lines = out.splitlines()
            
            # If the queue is empty or only contains headers, return zero consumption safely.
            if len(lines) <= 2:
                return 0, 0
                
            job_lines = lines[2:]
            total_threads = 0
            unique_samples = set()
            
            for line in job_lines:
                parts = line.split()
                if len(parts) < 8:
                    continue
                
                # Deduplicate sample IDs by assuming a standard job name syntax: '{SampleID}_{TaskIndex}_{TaskName}'
                job_name = parts[2]
                sample_id = job_name.split('_')[0]
                unique_samples.add(sample_id)
                
                # Extract the assigned slots (threads) located at the second to last index of standard SGE qstat layout.
                slots_str = parts[-2]
                slots = int(slots_str) if slots_str.isdigit() else 1
                total_threads += slots
                
            return len(unique_samples), total_threads
        except subprocess.CalledProcessError as e:
            # [MODIFIED] Upgraded from silent recovery (0,0) to a explicit error raise to catch cluster command drops.
            raise RuntimeError(f"[STRICT_ERROR] Failed to communicate with SGE cluster subsystem via qstat: {e.output}") from e
        except Exception as e:
            raise RuntimeError(f"[STRICT_ERROR] Unexpected structural parsing anomaly on qstat output: {str(e)}") from e

    def wait_for_resources(self, required_threads: int):
        """Blocks pipeline orchestration loop via polling until the required cluster slots are vacated."""
        # [MODIFIED] Added parameter safety guards against out-of-bounds task scaling requests.
        if required_threads is None or required_threads <= 0:
            raise ValueError("[STRICT_ERROR] Argument 'required_threads' must be a positive integer.")
        if required_threads > self.max_threads:
            raise ValueError(f"[STRICT_ERROR] Individual job requirement ({required_threads} threads) exceeds absolute global limit ({self.max_threads}).")

        while True:
            current_samples, current_threads = self._get_active_usage()
            
            sample_ok = current_samples < self.max_samples
            threads_ok = (current_threads + required_threads) <= self.max_threads
            
            if sample_ok and threads_ok:
                break
                
            print(f"[*] [Waiting] Resources Saturation: Samples ({current_samples}/{self.max_samples}), Threads ({current_threads}/{self.max_threads}). Polling sleep initiated (30s)...")
            time.sleep(30)

    def make_wrapper(self, cmd: str, wrapper_path: Path, work_dir: Path) -> Path:
        # [MODIFIED] Introduced execution safety checks to block creation of empty shell templates.
        if not cmd or not cmd.strip():
            raise ValueError("[STRICT_ERROR] Command execution line string 'cmd' cannot be empty.")
        if not wrapper_path or not isinstance(wrapper_path, Path):
            raise ValueError("[STRICT_ERROR] Destination 'wrapper_path' must be a valid Path object.")
        if not work_dir or not isinstance(work_dir, Path):
            raise ValueError("[STRICT_ERROR] Runtime execution tracking context 'work_dir' must be a valid Path object.")

        wrapper_path.parent.mkdir(parents=True, exist_ok=True)
        script = [
            "#!/bin/bash",
            "#$ -S /bin/bash",
            "set -e",
            f"cd {work_dir}", # [MODIFIED] Forced absolute path migration within shell context to keep relative file outputs bound correctly.
            cmd,
            f"echo OK > {wrapper_path.parent}/.done" # [MODIFIED] Unified filepath generation architecture using Path objects instead of os.path strings.
        ]
        wrapper_path.write_text("\n".join(script))
        wrapper_path.chmod(0o755)
        return wrapper_path

    # [MODIFIED] Eliminated default assignment for hold_jid. Callers must explicitly pass None if no dependency constraints exist.
    def qsub(self, job_id: str, script_path: Path, log_path: Path, threads: int, hold_jid: Any) -> str:
        # [MODIFIED] Added rigid pre-submission parameter integrity checks to bypass faulty command generation.
        if not job_id or not isinstance(job_id, str):
            raise ValueError("[STRICT_ERROR] Submission failed: 'job_id' must be a non-empty string.")
        if not script_path or not script_path.exists():
            raise FileNotFoundError(f"[STRICT_ERROR] Submission failed: Bash script asset not found at '{script_path}'.")
        if not log_path or not isinstance(log_path, Path):
            raise ValueError("[STRICT_ERROR] Submission failed: 'log_path' must be a valid Path directory object.")
        if threads is None or not isinstance(threads, int) or threads <= 0:
            raise ValueError("[STRICT_ERROR] Submission failed: Thread count allocation must be a positive integer.")

        log_path.mkdir(parents=True, exist_ok=True)

        # [MODIFIED] Switched to concrete logging nomenclature paths to bypass messy text replace side-effects.
        stdout_file = log_path / f"{job_id}.stdout"
        stderr_file = log_path / f"{job_id}.stderr"

        args = [
            "qsub", "-N", job_id, "-q", self.node,
            "-pe", "smp", str(threads),
            "-o", str(stdout_file), "-e", str(stderr_file),
            "-V", "-cwd"
        ]
        
        if hold_jid is not None:
            if not isinstance(hold_jid, str) or not hold_jid.strip():
                raise ValueError("[STRICT_ERROR] Validation failure: Inter-job dependency identifier 'hold_jid' must be a non-empty string.")
            args += ["-hold_jid", hold_jid]
        
        args.append(str(script_path))
        
        try:
            out = subprocess.check_output(args, text=True)
            return out.split()[2]
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"[STRICT_ERROR] SGE scheduler rejected job submission for '{job_id}': {e.output}") from e