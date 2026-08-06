"""
run.py - WPS 분석 CLI 진입점

분석 로직을 포함하지 않으며 argparse, 로깅 초기화,
모듈 호출 순서 정의만 담당한다.

실행 예시:
    python run.py \\
        --bam sample.bam \\
        --output-dir results/ \\
        --sample-name SAMPLE01 \\
        --ref-fasta hg38.fa \\
        --mappability-bw hg38.mappability.bw
"""

import argparse
import logging
import sys
from pathlib import Path

from utils import ensure_dir, get_logger, timer

# ── 패키지 수준 로거 (main에서 초기화) ────────────────────────────────────────
logger: logging.Logger = None  # type: ignore


# ── argparse ──────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    """CLI 인자 파서를 생성하여 반환한다."""
    parser = argparse.ArgumentParser(
        description="WPS (Window Protection Score) 계산 도구",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "예시:\n"
            "  python run.py --bam sample.bam --output-dir out/ --sample-name S01\n"
            "  python run.py --bam sample.bam --output-dir out/ --ref-fasta hg38.fa"
            " --mappability-bw hg38.map.bw"
        ),
    )

    # ── 필수 인자
    req = parser.add_argument_group("필수 인자")
    req.add_argument(
        "--bam", required=True, metavar="BAM",
        help="입력 BAM/CRAM 파일 경로 (인덱스 필요)"
    )
    req.add_argument(
        "--output-dir", required=True, metavar="DIR",
        help="결과 파일 출력 디렉토리 (없으면 자동 생성)"
    )

    # ── 선택 인자: 샘플
    samp = parser.add_argument_group("샘플 옵션")
    samp.add_argument(
        "--sample-name", default="sample", metavar="NAME",
        help="출력 파일명 prefix"
    )
    samp.add_argument(
        "--chroms", nargs="+", metavar="CHR",
        help="분석 대상 염색체 (기본: autosomal 전체, 예: chr1 chr2)"
    )

    # ── 선택 인자: Fragment 필터
    filt = parser.add_argument_group("Fragment 필터")
    filt.add_argument(
        "--min-frag-len", type=int, default=120, metavar="N",
        help="fragment 최소 길이 (bp)"
    )
    filt.add_argument(
        "--max-frag-len", type=int, default=180, metavar="N",
        help="fragment 최대 길이 (bp)"
    )
    filt.add_argument(
        "--min-mapq", type=int, default=20, metavar="Q",
        help="최소 매핑 품질 (MapQ)"
    )

    # ── 선택 인자: 정규화
    norm = parser.add_argument_group("정규화 옵션")
    norm.add_argument(
        "--ref-fasta", metavar="FA",
        help="GC content 보정용 참조 게놈 FASTA (미지정 시 GC 보정 생략)"
    )
    norm.add_argument(
        "--mappability-bw", metavar="BW",
        help="Mappability BigWig 파일 (미지정 시 mappability 마스킹 생략)"
    )
    norm.add_argument(
        "--mapq-threshold", type=float, default=0.9, metavar="F",
        help="Mappability 마스킹 임계값 (0.0–1.0)"
    )

    # ── 선택 인자: 출력
    out = parser.add_argument_group("출력 옵션")
    out.add_argument(
        "--no-bedgraph", action="store_true",
        help="bedGraph 파일 출력 생략"
    )
    out.add_argument(
        "--no-summary", action="store_true",
        help="통계 요약 TSV 출력 생략"
    )

    # ── 공통
    parser.add_argument(
        "--verbose", action="store_true",
        help="DEBUG 레벨 상세 로그 출력"
    )

    return parser


# ── 메인 오케스트레이터 ───────────────────────────────────────────────────────

def main(args: argparse.Namespace) -> None:
    """WPS 분석 파이프라인을 순서대로 실행한다.

    단계:
        1. 전처리: BAM 검증, 염색체 목록 확정
        2. 계산: short/long WPS 계산
        3. 정규화: GC 보정, mappability 마스킹
        4. 결과 저장: bedGraph, 통계 TSV 출력
    """
    global logger
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = get_logger("wps.run", log_level)

    logger.info("=== WPS 분석 시작 ===")
    logger.info(f"BAM: {args.bam}")
    logger.info(f"출력 디렉토리: {args.output_dir}")

    # 1. 전처리
    with timer("전처리", logger):
        from preprocess import run_preprocess
        prep = run_preprocess(
            bam_path=args.bam,
            min_frag_len=args.min_frag_len,
            max_frag_len=args.max_frag_len,
            min_mapq=args.min_mapq,
            chroms=args.chroms,
            logger=logger,
        )

    # 2. 계산
    with timer("WPS 계산", logger):
        from calculator import run_calculate
        wps_result = run_calculate(prep, logger=logger)

    # 3. 정규화
    with timer("정규화", logger):
        from normalize import run_normalize
        norm_result = run_normalize(
            wps_result=wps_result,
            ref_fasta=args.ref_fasta,
            mappability_bw=args.mappability_bw,
            mapq_threshold=args.mapq_threshold,
            logger=logger,
        )

    # 4. 결과 저장
    with timer("결과 저장", logger):
        from summary_results import run_summary
        run_summary(
            norm_result=norm_result,
            output_dir=args.output_dir,
            sample_name=args.sample_name,
            write_bedgraph=not args.no_bedgraph,
            write_summary_tsv=not args.no_summary,
            logger=logger,
        )

    logger.info("=== WPS 분석 완료 ===")


# ── 진입점 ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    try:
        main(args)
    except FileNotFoundError as e:
        print(f"[ERROR] 파일 없음: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"[ERROR] 파라미터 오류: {e}", file=sys.stderr)
        sys.exit(1)
    except ImportError as e:
        print(f"[ERROR] 패키지 누락: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] 사용자에 의해 중단됨", file=sys.stderr)
        sys.exit(130)