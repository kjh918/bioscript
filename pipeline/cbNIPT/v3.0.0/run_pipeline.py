import sys
import os
import argparse
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from src.bioscript.core.executor import SunGridExecutor
from src.bioscript.core.pipeline import Pipeline, Task


def get_parser():
    """CLI 옵션 파서를 생성하고 반환합니다. (매뉴얼 자동 생성기에서 참조)"""
    parser = argparse.ArgumentParser(description="cbNIPT Workflow Execution Script")
    
    # 필수 입출력 옵션
    parser.add_argument("-i", "--input", dest="raw_dir", required=True, help="Raw fastq.gz 파일이 있는 디렉토리 경로 (default: parsing fastq.gz)")
    parser.add_argument("-o", "--output", dest="work_dir", required=True, help="분석 결과가 저장될 디렉토리 경로")
    
    # 클러스터(SGE) 설정 옵션
    parser.add_argument("--sge-user", default="jhkim", help="SGE 실행 유저 이름 (기본값: jhkim)")
    parser.add_argument("--sge-node", default="all.q@ngsmaster", help="SGE 큐/노드 이름 (기본값: all.q@ngsmaster)")
    parser.add_argument("--max-samples", type=int, default=3, help="동시 실행할 최대 샘플 수 (기본값: 3)")
    parser.add_argument("--max-threads", type=int, default=24, help="최대 사용 스레드 수 (기본값: 24)")
    
    # 레퍼런스 설정 옵션 (기본값 세팅)
    parser.add_argument("--ref-fasta", default="/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa", help="Reference FASTA 경로")
    parser.add_argument("--known-snp", default="/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz", help="Known SNP VCF 경로")
    parser.add_argument("--known-indel1", default="/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz", help="Known Indel 1 VCF 경로")
    parser.add_argument("--known-indel2", default="/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", help="Known Indel 2 VCF 경로")


    return parser

if __name__ == "__main__":
    # --- [1] Argument Parser 설정 ---]
    parser = get_parser()
    args = parser.parse_args()
    
    # --- [2] 경로 및 변수 맵핑 ---
    RawFastqDir = args.raw_dir
    work_dir = Path(args.work_dir)
    scripts = Path(os.path.dirname(__file__)) / 'scripts'

    ReferenceFasta = args.ref_fasta
    KnownSnp = args.known_snp
    KnownIndel1 = args.known_indel1
    KnownIndel2 = args.known_indel2

    # --- [3] Pipeline Setting ---
    pipe = Pipeline(
        raw_dir = RawFastqDir,
        work_dir = work_dir,
        log_dir = work_dir / "logs",
        executor = SunGridExecutor(
            user=args.sge_user,
            node=args.sge_node,
            log_root= work_dir / "logs",
            max_samples=args.max_samples,  
            max_threads=args.max_threads  
        )
    )
        
    # --- [4] Tasks Setting ---
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
                    'KnownSnp': KnownSnp,
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
                runner_path= scripts / "run_gatk4_qc_alignment_summary_using_singularity.py",
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
                runner_path= scripts / "run_gatk4_qc_insert_size_using_singularity.py",
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
                runner_path= scripts / "runs_gatk4_qc_wgs_metrics_using_singularity.py",
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
                log_path = work_dir / sid / "logs" / "12_mosdepth_singularity_coverage_100kb", 
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