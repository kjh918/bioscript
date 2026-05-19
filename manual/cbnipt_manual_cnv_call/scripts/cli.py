import argparse
import sys
import os
from pipeline import CnvPipeline
from utils import log

def main():
    parser = argparse.ArgumentParser(description="cbNIPT & Single-cell WGA CNV Analyzer (Single BAM)")

    # 1. 필수 입력 정보
    group_input = parser.add_argument_group('Input/Output')
    group_input.add_argument('--SeqID', required=True, help='분석 고유 ID (파일명에 사용)')
    group_input.add_argument('--BamPath', required=True, help='분석할 단일 BAM 파일 경로')
    group_input.add_argument('--ReferenceFasta', required=True, help='참조 유전체 FASTA 파일 경로')
    group_input.add_argument('--OutDir', required=True, help='결과 저장 디렉토리')

    # 2. Analysis Parameters
    group_params = parser.add_argument_group('Analysis Parameters')
    group_params.add_argument('--BinSize', type=int, default=100000, help='Bin 크기 (100kb)')
    group_params.add_argument('--MappabilityBW', default='/storage/home/jhkim/Projects/cbNIPT/260423-GCX-cbNIPT-ManualMethod/Resources/reference/hg38.100mer.bw', help='Mappability BigWig 파일')
    group_params.add_argument('--MinMapQ', type=int, default=20, help='최소 Mapping Quality')
    group_params.add_argument('--LowessFrac', type=float, default=0.2, help='GC 보정 LOWESS 비율')
    group_params.add_argument('--BaselinePloidy', type=int, default=2, help='기준 Ploidy')
    # [추가] 스무딩 윈도우 (데이터 노이즈 제거용)
    group_params.add_argument('--SmoothWindow', type=int, default=5, help='Rolling median window size')

    # 3. Filtering & Normalization
    group_filter = parser.add_argument_group('Filtering & Normalization')
    group_filter.add_argument('--MinMappability', type=float, default=0.9, help='최소 Mappability 점수')
    group_filter.add_argument('--MinDepth', type=float, default=2, help='최소 Depth')
    group_filter.add_argument('--MinCoverage', type=float, default=0.5, help='최소 Coverage 비율')
    group_filter.add_argument('--IncludeSexChrom', type=bool, default=True, help='성염색체 포함 여부')
    group_filter.add_argument('--SegPenalty', type=float, default=10.0, help='Segmentation 민감도')
    
    # [추가] 3중 필터링을 위한 상세 설정 인자들
    group_filter.add_argument('--MinGC', type=float, default=0.3, help='최소 GC 함량 필터 (기본 0.3)')
    group_filter.add_argument('--MaxGC', type=float, default=0.7, help='최대 GC 함량 필터 (기본 0.7)')
    group_filter.add_argument('--MaskLowerQ', type=float, default=0.01, help='하위 Coverage 필터 (기본 0.01)')
    group_filter.add_argument('--MaskUpperQ', type=float, default=0.99, help='상위 Coverage 필터 (기본 0.99)')
    group_filter.add_argument('--Threads', type=float, default=4, help='상위 Coverage 필터 (기본 0.99)')

    args = parser.parse_args()

    # 분석 시작 로그
    log("="*50)
    log(f" Starting scDNA/cbNIPT Analysis pipeline")
    log(f" Run ID     : {args.SeqID}")
    log(f" BAM Source : {args.BamPath}")
    log(f" Bin Size   : {args.BinSize / 1000} kb")
    log("="*50)

    try:
        # 4. 파이프라인 객체 생성 및 가동
        pipeline = CnvPipeline(args)
        pipeline.run()
        
        log(f"--- [SUCCESS] {args.SeqID} Analysis finished successfully ---")

    except Exception as e:
        log(f"--- [ERROR] Analysis failed: {str(e)} ---")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()