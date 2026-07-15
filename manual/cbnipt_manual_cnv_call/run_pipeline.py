import sys
import os
import argparse
from pathlib import Path

sys.path.append(os.path.dirname(__file__))

from src.bioscript.core.executor import SunGridExecutor
from src.bioscript.core.pipeline import Pipeline, Task

def get_parser():
    """CLI 옵션 파서를 생성하고 반환합니다."""
    parser = argparse.ArgumentParser(description="cbNIPT Workflow Execution Script")
    parser.add_argument("-i", "--input", dest="raw_dir", required=True, help="Raw fastq/bam 디렉토리 경로")
    parser.add_argument("-o", "--output", dest="work_dir", required=True, help="결과 저장 디렉토리")
    parser.add_argument("--sge-user", default="jhkim", help="SGE 유저 이름")
    parser.add_argument("--sge-node", default="all.q@ngsmaster", help="SGE 큐/노드")
    parser.add_argument("--max-samples", type=int, default=3, help="동시 실행 최대 샘플 수")
    parser.add_argument("--max-threads", type=int, default=24, help="최대 스레드 수")
    return parser

if __name__ == "__main__":
    # --- [1] Argument Parser 설정 ---
    parser = get_parser()
    args = parser.parse_args()
    
    # --- [2] 경로 및 변수 맵핑 ---
    RawBamDir = args.raw_dir
    work_dir = Path(args.work_dir)
    
    # [MODIFIED] script(생성된 래퍼)와 components(원본 YAML) 경로를 명확히 분리
    scripts_dir = Path(os.path.dirname(__file__)) / 'scripts'

    # --- [3] Pipeline Setting ---
    pipe = Pipeline(
        raw_dir = RawBamDir,
        work_dir = work_dir,
        log_dir = work_dir / "logs",
        executor = SunGridExecutor(
            user=args.sge_user,
            node=args.sge_node,
            log_root=work_dir / "logs",
            max_samples=args.max_samples,  
            max_threads=args.max_threads  
        ),
        suffix='bam'
    )
    
    # --- [4] Global Task Setting (샘플과 무관하게 1회만 실행) ---
    global_task = Task(
        name="Global_Reference_Binning",
        runner_path = scripts_dir / "run_cbnipt_makebins_generate_reference_genome_bins.py",
        log_path = work_dir / "logs" / "00.global_binning", 
        spec = {
            'script_path': scripts_dir / "cli.py",
            'ReferenceFasta': "/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa",
            'MappabilityBW': "/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/hg38.100mer.bw",
            'OutBinFile': work_dir / "Binning",
            'BinSize': 100000,
            'MinMappability': 0.9,
            'MinGC': 0.3,
            'MaxGC': 0.7,
            "Threads": 1
        }
    )
    pipe.add_tasks([global_task])
    ## --- [5] Sample Tasks Setting ---
    
    for sid in pipe.samples:
        sample_tasks = [
            Task(
                name="CNVCall",
                runner_path = scripts_dir / "run_cbnipt_callcnv_calculate_cnv_and_variant_profile.py", # 적절한 스크립트로 변경 필요
                log_path = work_dir / sid / "logs" / "02.cnv_call", 
                spec = {
                    "script_path":scripts_dir / "cli.py",
                    'SeqID': sid,
                    
                    ##  check options
                    'MinCoverage':0,

                    'BamPath': f"{RawBamDir}/{sid}.bam",
                    'VcfFile':'/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/KOVA_v7/kova_sites_vcf/KOVA_v7_merged.vcf.gz',
                    'AnnotatedBins': work_dir / "Binning" / "hg38.fa.bin_100.0K.bed.gz",
                    'OutDir': work_dir / sid,
                    
                    "Threads": 4
                }
            )
        ]
        # 개별 샘플 Task들을 파이프라인에 등록
        pipe.add_tasks(sample_tasks)

    # --- [6] Execute ---
    # Global Task가 먼저 qsub 되고, 그게 끝나야(-hold_jid) Sample Task들이 병렬로 실행됨
    pipe.run(run_integrity=True)