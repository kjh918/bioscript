import argparse
import sys
import os
import pysam
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))  # 프로젝트 루트 디렉토리 추가

from src.cbnipt_cnv_caller.pipeline import CnvPipeline
from src.cbnipt_cnv_caller.binning import generate_bins, annotate_bin_metadata, get_chromosomes, apply_final_filters

def run_make_bins(args):
    """
    [독립 실행 모듈]
    레퍼런스 FASTA를 기반으로 Bin을 쪼개고, GC/Mappability 필터링을 거친 후
    고속 쿼리가 가능한 압축/인덱싱된 BED 파일(.bed.gz 및 .tbi)을 생성합니다.
    """
    print("="*50)
    print(f" [Step 0] Generating Pre-computed Reference Bins (BED.GZ)")
    print(f" Bin Size : {args.BinSize / 1000} kb")
    print(f" Filters  : GC({args.MinGC}~{args.MaxGC}), MapQ>={args.MinMappability}")
    print("="*50)
    
    # 1. Binning & Annotation
    chroms = get_chromosomes(fasta_path=args.ReferenceFasta, include_sex=args.IncludeSexChrom)
    bins = generate_bins(fasta_path=args.ReferenceFasta, bin_size=args.BinSize, chromosomes=chroms)
    bins_annotated = annotate_bin_metadata(bins, fasta_path=args.ReferenceFasta, mappability_bw=args.MappabilityBW)
    
    # 2. Filtering (binning.py의 apply_final_filters 호출)
    filtered_bins = apply_final_filters(bins_annotated, args)
    
    # 3. BED 표준 포맷 준수 및 Tabix 인덱싱
    cols = ['chrom', 'start', 'end'] + [c for c in filtered_bins.columns if c not in ['chrom', 'start', 'end']]
    filtered_bins = filtered_bins[cols]
    
    # 헤더 '#' 처리 (IGV 및 tabix 호환성)
    filtered_bins.rename(columns={'chrom': '#chrom'}, inplace=True)
    
    temp_out = f'{args.OutBinFile}/{os.path.basename(args.ReferenceFasta)}.bin_{str(round(args.BinSize/1000,0))}K.bed'
    os.makedirs(os.path.dirname(os.path.abspath(temp_out)), exist_ok=True)
    filtered_bins.to_csv(temp_out, sep="\t", index=False)
    
    # bgzip 압축 및 tabix 인덱싱 실행
    pysam.tabix_index(temp_out, preset='bed', force=True)
    
    if os.path.exists(temp_out):
        os.remove(temp_out)
        
    final_output = f"{temp_out}.gz"
    print(f"--- [SUCCESS] Indexed BED saved to: {final_output} ---")
    print(f"--- [SUCCESS] Tabix Index saved to: {final_output}.tbi ---")


def run_cnv_pipeline(args):
    """
    [메인 분석 모듈]
    미리 생성된 Bin 파일을 입력으로 받아 단일 샘플의 BAM 분석 파이프라인을 가동합니다.
    """
    print("="*50)
    print(f" Starting cbNIPT Analysis pipeline")
    print(f" Run ID      : {args.SeqID}")
    print(f" BAM Source  : {args.BamPath}")
    print(f" Target Bins : {args.AnnotatedBins}")
    print("="*50)

    try:
        pipeline = CnvPipeline(args)
        pipeline.run()
        print(f"--- [SUCCESS] {args.SeqID} Analysis finished successfully ---")

    except Exception as e:
        print(f"--- [ERROR] Analysis failed: {str(e)} ---")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="cbNIPT & Single-cell WGA CNV Analyzer")
    subparsers = parser.add_subparsers(dest="command", required=True, help="실행할 파이프라인 모드 선택")

    # ---------------------------------------------------------
    # 모드 1: make-bins (레퍼런스 고정 빈 생성기)
    # ---------------------------------------------------------
    parser_bins = subparsers.add_parser('make-bins', help='Generate pre-computed and annotated reference bins')
    
    parser_bins.add_argument('--ReferenceFasta', required=True, help='참조 유전체 FASTA 파일 경로')
    parser_bins.add_argument('--MappabilityBW', required=True, help='Mappability BigWig 파일')
    parser_bins.add_argument('--OutBinFile', required=True, help='생성된 빈 정보를 저장할 경로 (예: hg38_100kb_bins.bed.gz)')
    
    parser_bins.add_argument('--BinSize', type=int, default=100000, help='Bin 크기 (100kb)')
    parser_bins.add_argument('--MinMappability', type=float, default=0.9, help='최소 Mappability 점수')
    parser_bins.add_argument('--MinGC', type=float, default=0.3, help='최소 GC 함량 필터')
    parser_bins.add_argument('--MaxGC', type=float, default=0.7, help='최대 GC 함량 필터')
    
    # 문자열 입력을 boolean으로 처리하기 위한 트릭
    parser_bins.add_argument('--IncludeSexChrom', action='store_true', default=True, help='성염색체 포함 여부')

    # ---------------------------------------------------------
    # 모드 2: call-cnv (단일 샘플 CNV 파이프라인)
    # ---------------------------------------------------------
    parser_cnv = subparsers.add_parser('call-cnv', help='Run CNV pipeline on a single BAM using pre-computed bins')
    
    parser_cnv.add_argument('--SeqID', required=True, help='분석 고유 ID')
    parser_cnv.add_argument("--VcfFile", type=str, required=True, help="Path to the merged VCF file (.vcf.gz)")
    parser_cnv.add_argument('--BamPath', required=True, help='분석할 단일 BAM 파일 경로')
    parser_cnv.add_argument('--AnnotatedBins', required=True, help='make-bins로 생성한 사전 계산 빈 파일(.bed.gz)')
    parser_cnv.add_argument('--OutDir', required=True, help='결과 저장 디렉토리')
    
    parser_cnv.add_argument('--MinMapQ', type=int, default=20, help='최소 Mapping Quality')
    parser_cnv.add_argument('--LowessFrac', type=float, default=0.2, help='GC 보정 LOWESS 비율')
    parser_cnv.add_argument('--BaselinePloidy', type=int, default=2, help='기준 Ploidy')
    parser_cnv.add_argument('--SmoothWindow', type=int, default=5, help='Rolling median window size')
    parser_cnv.add_argument('--MinDepth', type=float, default=2, help='최소 Depth')
    parser_cnv.add_argument('--MinCoverage', type=float, default=0.5, help='최소 Coverage 비율')
    parser_cnv.add_argument('--MaskLowerQ', type=float, default=0.01, help='하위 Coverage 필터')
    parser_cnv.add_argument('--MaskUpperQ', type=float, default=0.99, help='상위 Coverage 필터')
    parser_cnv.add_argument('--SegPenalty', type=float, default=10.0, help='Segmentation 민감도')
    parser_cnv.add_argument('--Threads', type=int, default=4, help='사용할 스레드 수')

    args = parser.parse_args()

    if args.command == 'make-bins':
        run_make_bins(args)
    elif args.command == 'call-cnv':
        run_cnv_pipeline(args)

if __name__ == "__main__":
    main()