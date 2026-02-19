
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from src.bioscript.core.executor import SunGridExecutor
from src.bioscript.core.pipeline import Pipeline, Task
from pathlib import Path
from typing import Dict, Any, List

if __name__ == "__main__":

    ## Workflow Setting ## 
    work_dir = Path("/storage/home/jhkim/Projects/cbNIPT/260122-GCX-cbNIPT_Workflow/Results")
    scripts = Path(os.path.dirname(__file__)) / 'scripts'

    ## Resources Setting ##
    RawFastqDir= '/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Resources/Rawdata/PicoPLEXGold/251103'
    ReferenceFasta='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa'
    KnownSnp='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz'
    KnownIndel1='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz'
    KnownIndel2='/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
        
    ## Pipeline Setting ##
    pipe = Pipeline(
        raw_dir = RawFastqDir,
        work_dir = work_dir,
        log_dir = work_dir / "logs",
        executor = SunGridExecutor(
            user='jhkim',
            node="all.q@ngsmaster",
            log_root= work_dir / "logs",
            max_samples=3,  
            max_threads=24  
        )
    )
        
    ## Tasks Setting ##
    for sid in pipe.samples:
        tasks = [
            Task(
                name="fastqc",
                runner_path= scripts / "run_fastqc_singularity_pe_qc_extract.py",
                log_path = work_dir / sid / "logs" / "00_fastqc", 
                spec = {
                    'SeqID': sid,
                    'RawFastqDir': RawFastqDir,
                    'qcResDir': work_dir / sid / "fastq_raw",
                    "Threads": 4,
                }
            ),
            Task(
                name="fastp",
                runner_path= scripts / "run_fastp_picoplexgold_pe_trim_using_singularity.py",
                log_path = work_dir / sid / "logs" / "01_fastp", 
                spec = {
                    'SeqID': sid,
                    'RawFastqDir': RawFastqDir,
                    'TrimFastqDir': work_dir / sid / "fastq_trimmed",
                    'qcResDir': work_dir / sid / "qc",
                    "Threads": 4,
                }
            ),
            Task(
                name="bwa_picard_align_merge",
                runner_path= scripts / "run_bwa_picard_pe_align_merge.py",
                log_path = work_dir / sid / "logs" / "02_bwa_picard_align_merge", 
                spec = {
                    'SeqID': sid,
                    'TrimFastqDir': work_dir / sid / "fastq_trimmed",
                    'BamDir': work_dir / sid / "bam",
                    'TmpDir': work_dir / sid / "tmp",
                    'ReferenceFasta': ReferenceFasta,
                    'ReadGroupID': sid,
                    'ReadGroupPlatform': 'ILLUMINA',
                    'ReadGroupLibrary': 'PicoPLEXGold',
                    'ReadGroupCenter': 'GCX',
                    "Threads": 8,
                }
            ),
            Task(
                name="samtools_bam_sort_index",
                runner_path= scripts / "run_samtools_bam_sort_index.py",
                log_path = work_dir / sid / "logs" / "03_samtools_bam_sort_index", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    "Threads": 2,
                }
            ),
            Task(
                name="gatk4_singularity_dedup",
                runner_path= scripts / "run_gatk4_dedup_using_singularity.py",
                log_path = work_dir / sid / "logs" / "04_gatk4_singularity_dedup", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk3_singularity_indel_realign",
                runner_path= scripts / "run_gatk3_indel_realign_using_singularity.py",
                log_path = work_dir / sid / "logs" / "05_gatk3_singularity_indel_realign", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    'KnownIndel1': KnownIndel1,
                    'KnownIndel2': KnownIndel2,
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk4_singularity_bqsr",
                runner_path= scripts / "run_gatk4_bqsr_using_singularity.py",
                log_path = work_dir / sid / "logs" / "06_gatk4_singularity_bqsr", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    'KnownIndel1': KnownIndel1,
                    'KnownIndel2': KnownIndel2,
                    "Threads": 8,
                }
            ),
            Task(
                name="samtools_filter_index_recal_bam_filter",
                runner_path= scripts / "run_samtools_filter_bam.py",
                log_path = work_dir / sid / "logs" / "07_samtools_filter_index_recal_bam_filter", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk4_singularity_qc_artifacts",
                runner_path= scripts / "run_gatk4_qc_artifacts_using_singularity.py",
                log_path = work_dir / sid / "logs" / "08_gatk4_singularity_qc_artifacts", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk4_singularity_qc_alignment_summary",
                runner_path= scripts / "09_gatk4_singularity_qc_alignment_summary.py",
                log_path = work_dir / sid / "logs" / "09_gatk4_singularity_qc_alignment_summary", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk4_singularity_qc_insert_size",
                runner_path= scripts / "run_gatk4_qc_insert_size_using_singularity.py*",
                log_path = work_dir / sid / "logs" / "10_gatk4_singularity_qc_insert_size", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
            Task(
                name="gatk4_singularity_qc_wgs_metrics",
                runner_path= scripts / "run_gatk4_qc_wgs_metrics_using_singularity.py*",
                log_path = work_dir / sid / "logs" / "11_gatk4_singularity_qc_wgs_metrics", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
            Task(
                name="mosdepth_singularity_coverage_100kb",
                runner_path= scripts / "run_mosdepth_singularity_coverage.py",
                log_path = work_dir / sid / "logs" / "11_mosdepth_singularity_coverage_100kb", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid /"bam",
                    'qcResDir': work_dir / sid / "qc",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
            Task(
                name="copykit",
                runner_path= scripts / "run_copykit_cbnipt_copynumber_analysis.py",
                log_path = work_dir / sid / "logs" / "13_copykit_cbnipt_copynumber_analysis.py", 
                spec = {
                    'SeqID': sid,
                    'NGS_DataBaseDir': work_dir / sid / "bam",
                    'ResultBaseDir': work_dir / sid / "copykit",
                    'ReferenceFasta': ReferenceFasta,
                    "Threads": 8,
                }
            ),
        ]
        
        pipe.add_tasks(tasks)
        
    pipe.run()

