"""
run.py — nipt_fragmentomics CLI entry point.

사용법
------
# 기본 (CNV + FF만)
python -m run \\
    --bam sample.bam \\
    --out-dir ./results/SID001

# WPS 포함 (L-WPS, long fragment)
python -m run \\
    --bam sample.bam \\
    --out-dir ./results/SID001 \\
    --marker-bed cCREs.bed \\
    --wps-mode L \\
    --wps-frag long

# 전체 옵션
python -m run \\
    --bam      sample.bam \\
    --out-dir  ./results/SID001 \\
    --bin-bed  bins_100kb.bed.gz \\
    --marker-bed cCREs.bed \\
    --fasta    hg38.fa \\
    --bw       hg38.100mer.bw \\
    --vcf      snp.vcf.gz \\
    --wps-mode L S \\
    --wps-frag long short \\
    --wps-extend 2000 \\
    --bin-size 100000 \\
    --jobs 8 \\
    --resume \\
    --sample-id SID001
"""

from __future__ import annotations

import argparse
import logging
import os
import sys


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run.py",
        description="NIPT fragmentomics 통합 파이프라인",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 필수
    p.add_argument("--bam",      required=True, help="인덱스된 BAM 파일")
    p.add_argument("--out-dir",  required=True, dest="out_dir",
                   help="결과 저장 디렉토리")

    # 참조 파일
    ref = p.add_argument_group("참조 파일 (선택)")
    ref.add_argument("--bin-bed",    dest="bin_bed",    default=None,
                     help="bin 영역 BED — 없으면 --bin-size 로 자동 생성")
    ref.add_argument("--marker-bed", dest="marker_bed", default=None,
                     help="Marker BED (WPS profile용, chrom start end marker_id marker_type)")
    ref.add_argument("--fasta",  default=None, help="GC 계산용 참조 FASTA (.fai 필요)")
    ref.add_argument("--bw",     default=None, help="Mappability bigWig")
    ref.add_argument("--vcf",    default=None,
                     help="Population SNP VCF (BAF 계산용, .tbi 인덱스 필요)")

    # bin 파라미터
    bp = p.add_argument_group("Bin 파라미터")
    bp.add_argument("--bin-size",  dest="bin_size", type=int, default=100_000,
                    help="CNV용 자동 bin 크기 (bp), bin-bed 없을 때 사용")
    bp.add_argument("--min-mapq",  dest="min_mapq",  type=int, default=20)
    bp.add_argument("--min-baseq", dest="min_baseq", type=int, default=20)
    bp.add_argument("--min-mappability", dest="min_mappability",
                    type=float, default=0.75)

    # WPS 파라미터
    wp = p.add_argument_group("WPS 파라미터 (미지정 시 WPS 생략)")
    wp.add_argument("--wps-mode",   dest="wps_modes", nargs="*",
                    choices=["L", "S"], default=None,
                    help="WPS 모드: L=뉴클레오솜(k=120bp) S=TF footprint(k=16bp). "
                         "예: --wps-mode L S")
    wp.add_argument("--wps-frag",   dest="wps_frags", nargs="*",
                    choices=["short", "long"], default=None,
                    help="fragment 필터 (기본: long). 예: --wps-frag long short")
    wp.add_argument("--wps-extend", dest="wps_extend", type=int, default=2000,
                    help="marker profile 추출 시 중심 ±extend bp")

    # CNV 파라미터
    cp = p.add_argument_group("CNV 파라미터")
    cp.add_argument("--zscore-gain", dest="zscore_gain", type=float, default=3.0)
    cp.add_argument("--zscore-loss", dest="zscore_loss", type=float, default=-3.0)

    # BAF 파라미터
    baf = p.add_argument_group("BAF 파라미터 (--vcf 지정 시 적용)")
    baf.add_argument("--baf-af-min",    type=float, default=0.2,  dest="baf_af_min")
    baf.add_argument("--baf-af-max",    type=float, default=0.8,  dest="baf_af_max")
    baf.add_argument("--baf-min-depth", type=int,   default=5,    dest="baf_min_depth",
                     help="site 당 최소 depth")

    # 실행 옵션
    ex = p.add_argument_group("실행 옵션")
    ex.add_argument("--jobs",       type=int,  default=4,
                    help="염색체 병렬 처리 수")
    ex.add_argument("--resume",     action="store_true",
                    help="중간 파일 있으면 해당 Step 건너뜀")
    ex.add_argument("--no-viz",     action="store_true",
                    help="Plotly 시각화 파일 생성 생략")
    ex.add_argument("--sample-id",  dest="sample_id", default="",
                    help="샘플 ID (로그·그래프 제목에 표시)")
    ex.add_argument("--log-level",  dest="log_level", default="INFO",
                    choices=["DEBUG","INFO","WARNING","ERROR"])

    return p


def _setup_logging(out_dir: str, level: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    fmt = "%(asctime)s [%(levelname)s] %(name)s — %(message)s"
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format=fmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(os.path.join(out_dir, "pipeline.log"), mode="a"),
        ],
    )


def main() -> None:
    parser = _build_parser()
    args   = parser.parse_args()

    # 파일 존재 확인
    if not os.path.exists(args.bam):
        parser.error(f"BAM 파일을 찾을 수 없습니다: {args.bam}")
    for opt, path in [("--bin-bed",    args.bin_bed),
                      ("--marker-bed", args.marker_bed),
                      ("--fasta",      args.fasta),
                      ("--bw",         args.bw),
                      ("--vcf",        args.vcf)]:
        if path and not os.path.exists(path):
            parser.error(f"{opt} 파일을 찾을 수 없습니다: {path}")

    _setup_logging(args.out_dir, args.log_level)

    from nipt_fragmentomics.pipeline import run
    run(
        bam_path         = args.bam,
        out_dir          = args.out_dir,
        bin_bed          = args.bin_bed,
        marker_bed       = args.marker_bed,
        fasta_path       = args.fasta,
        bw_path          = args.bw,
        vcf_path         = args.vcf,
        bin_size         = args.bin_size,
        min_mapq         = args.min_mapq,
        min_baseq        = args.min_baseq,
        min_mappability  = args.min_mappability,
        wps_modes        = args.wps_modes,
        wps_frags        = args.wps_frags,
        wps_extend       = args.wps_extend,
        zscore_gain      = args.zscore_gain,
        zscore_loss      = args.zscore_loss,
        baf_af_min       = args.baf_af_min,
        baf_af_max       = args.baf_af_max,
        baf_min_depth    = args.baf_min_depth,
        n_jobs           = args.jobs,
        resume           = args.resume,
        make_viz         = not args.no_viz,
        sample_id        = args.sample_id,
    )


if __name__ == "__main__":
    main()