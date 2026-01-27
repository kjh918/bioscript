
from src.executor import SunGridExecutor
from src.pipeline import Pipeline, Task
from pathlib import Path
import os
from typing import Dict, Any, List

if __name__ == "__main__":

	work_dir = Path("/storage/home/jhkim/pipelines/cbNIPT/v0.0.3/work")
	scripts = Path("/storage/home/jhkim/pipelines/cbNIPT/v0.0.3/scripts")

	pipe = Pipeline(
		raw_dir = '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Resources/Rawdata/PicoPLEXGold/251103',
		work_dir = work_dir,
		log_dir = work_dir / "logs",
		executor = SunGridExecutor(
			user='jhkim',
			node="all.q@ngsnode1",
			log_root= work_dir / "logs",
			max_samples=3,	# 동시에 최대 3개 샘플만 진행
			max_threads=24	# 총 사용 스레드가 24개를 넘지 않도록 제어
		)
	)

	for sid in pipe.samples:
		tasks = [
			Task(
				name="fastqc",
				runner_path= scripts / "run_fastqc_singularity_pe_qc_extract.py",
				log_path = work_dir / sid / "fastq_raw" / "logs", 
				spec = {
					'SeqID': sid,
					'RawFastqDir': '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Resources/Rawdata/PicoPLEXGold/251103',
					'qcResDir': work_dir / sid / "fastq_raw",
					"threads": 8,
				}
			),
			#Task(
			#	name="fastp",
			#	runner_path= scripts / "run_fastqc_0_12_1_singularity_pe_qc_extract.py",
			#	log_path = workf_dir / sid / "01_fastqc" / "logs", 
			#	spec = {
			#		'SeqID': sid,
			#		'RawFastqDir': /storage/home/jhkim/pipelines/cbNIPT/v0.0.3/work,
			#		'TrimFastqDir': work_dir / "trimmed",
			#		'qcResDir': work_dir / "qc",
			#		"threads": 8,
			#	}
			#),
		]
		
		pipe.add_tasks(tasks)
	
	pipe.run()

