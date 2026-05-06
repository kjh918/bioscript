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
    parser.add_argument("--sge-node", default="all.q@ngsnode1", help="SGE 큐/노드 이름 (기본값: all.q@ngsmaster)")
    parser.add_argument("--max-samples", type=int, default=3, help="동시 실행할 최대 샘플 수 (기본값: 3)")
    parser.add_argument("--max-threads", type=int, default=24, help="최대 사용 스레드 수 (기본값: 24)")
    
    # 레퍼런스 설정 옵션 (기본값 세팅)
    #parser.add_argument("--ref-fasta", default="/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa", help="Reference FASTA 경로")
    #parser.add_argument("--known-snp", default="/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz", help="Known SNP VCF 경로")
    #parser.add_argument("--known-indel1", default="/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz", help="Known Indel 1 VCF 경로")
    #parser.add_argument("--known-indel2", default="/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", help="Known Indel 2 VCF 경로")


    return parser

if __name__ == "__main__":
    # --- [1] Argument Parser 설정 ---]
    parser = get_parser()
    args = parser.parse_args()
    
    # --- [2] 경로 및 변수 맵핑 ---
    RawFastqDir = args.raw_dir
    work_dir = Path(args.work_dir)
    scripts = Path(os.path.dirname(__file__)) / 'scripts'

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
    for sid in pipe.samples[:1]:
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
                runner_path= scripts / "run_fastp_pe_trim_using_singularity.py",
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
                name="star",
                runner_path= scripts / "run_star_alignment_rna_seq_mapping.py",
                log_path = work_dir / sid / "logs" / "02_star", 
                spec = {
                    'SeqID': sid,
                    'FastqDir': work_dir / sid / "fastq_trimmed",
                    'BamDir': work_dir / sid / "bam",
                    #'StarIndex': args.star_index,
                    'StarIndex': "/storage/references_and_index/hg38/star-rsem-index/PE-150",
                    "Threads": 20,
                }
            ),
            Task(
                name="rsem",
                runner_path= scripts / "run_rsem_calculate_expression_quanti.py",
                log_path = work_dir / sid / "logs" / "03_rsem", 
                spec = {
                    'SeqID': sid,
                    'TranscriptomeBam': work_dir / sid / "bam",
                    'RsemIndex': "/storage/references_and_index/hg38/star-rsem-index/PE-150/hg38",
                    'OutputDir': work_dir / sid / "quanti" / "rsem", 
                    "Threads": 20,
                }
            ),
            Task(
                name="salmon",
                runner_path= scripts / "run_salmon_transcript_quanti.py",
                log_path = work_dir / sid / "logs" / "04_salmon", 
                spec = {
                    'SeqID': sid,
                    'BamDir': work_dir / sid / "bam",
                    'SalmonIndex': "/storage/references_and_index/hg38/star-rsem-index/PE-150/hg38.transcripts.fa",
                    'OutputDir': work_dir / sid / "quanti" / "salmon", 
                    "Threads": 20,
                }
            ),





            #Task(
            #    name="bbsplit",
            #    runner_path= scripts / "run_bwa_picard_pe_align_merge.py",
            #    log_path = work_dir / sid / "logs" / "02_bwa_picard_align_merge", 
            #    spec = {
            #        'SeqID': sid,
            #        'TrimFastqDir': work_dir / sid / "fastq_trimmed",
            #        'BamDir': work_dir / sid / "bam",
            #        'TmpDir': work_dir / sid / "tmp",
            #        'ReferenceFasta': ReferenceFasta,
            #        'ReadGroupID': sid,
            #        'ReadGroupPlatform': 'ILLUMINA',
            #        'ReadGroupLibrary': 'PicoPLEXGold',
            #        'ReadGroupCenter': 'GCX',
            #        "Threads": 8,
            #    }
            #),
            #Task(
            #    name="sortmerna",
            #    runner_path= scripts / "run_bwa_picard_pe_align_merge.py",
            #    log_path = work_dir / sid / "logs" / "02_bwa_picard_align_merge", 
            #    spec = {
            #        'SeqID': sid,
            #        'TrimFastqDir': work_dir / sid / "fastq_trimmed",
            #        'BamDir': work_dir / sid / "bam",
            #        'TmpDir': work_dir / sid / "tmp",
            #        'ReferenceFasta': ReferenceFasta,
            #        'ReadGroupID': sid,
            #        'ReadGroupPlatform': 'ILLUMINA',
            #        'ReadGroupLibrary': 'PicoPLEXGold',
            #        'ReadGroupCenter': 'GCX',
            #        "Threads": 8,
            #    }
            #),

        ]
        
        pipe.add_tasks(tasks)
        
    pipe.run()