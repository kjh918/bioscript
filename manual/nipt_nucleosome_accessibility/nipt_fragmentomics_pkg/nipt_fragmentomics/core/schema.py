"""
schema.py
=========
Fragment 단위 점수를 담는 경량 dataclass.

BAM 스캔 시 각 fragment 마다 생성되며,
bin accumulator 에 집계된 뒤 즉시 소멸합니다.
DataFrame 변환 없이 직접 bin 누적기에 전달하므로
메모리·속도 오버헤드가 최소입니다.
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(slots=True)
class FragmentScore:
    """
    단일 cfDNA fragment 의 측정 결과.

    Attributes
    ----------
    qname       : read name (fragment 식별자)
    chrom       : 염색체
    frag_start  : 0-based fragment 시작 좌표 (inclusive)
    frag_end    : 0-based fragment 끝 좌표 (exclusive)
    frag_len    : fragment 길이 (bp)
    is_short    : True → short (≤150 bp), False → long (≥151 bp)
    midpoint    : (frag_start + frag_end) // 2  — bin 귀속 기준
    mapq        : MAPQ 값
    """
    qname:      str
    chrom:      str
    frag_start: int
    frag_end:   int
    frag_len:   int
    is_short:   bool
    midpoint:   int
    mapq:       int

    @classmethod
    def from_read(cls, read, max_frag_len: int = 1000) -> "FragmentScore | None":
        """
        pysam AlignedSegment 로부터 FragmentScore 를 생성합니다.

        Parameters
        ----------
        max_frag_len : 허용 최대 fragment 길이 (bp).
                       비정상 template_length (수백만 bp) 필터링.
                       cfDNA 분석에서 1000bp 이상은 일반적으로 비정상.
        """
        if read.is_paired and read.template_length != 0:
            tlen = abs(read.template_length)
            # 비정상 template_length 필터: read 길이보다 작거나 너무 큰 경우
            # next_reference_start 기반으로 실제 insert size 재계산
            if tlen > max_frag_len or tlen < (read.query_length or 0):
                # template_length 를 신뢰할 수 없으면 reference_end 사용
                frag_start = read.reference_start
                frag_end   = read.reference_end
            else:
                frag_start = read.reference_start
                frag_end   = frag_start + tlen
        else:
            frag_start = read.reference_start
            frag_end   = read.reference_end

        frag_len = frag_end - frag_start
        if frag_len <= 0 or frag_len > max_frag_len:
            return None

        return cls(
            qname      = read.query_name,
            chrom      = read.reference_name,
            frag_start = frag_start,
            frag_end   = frag_end,
            frag_len   = frag_len,
            is_short   = (frag_len <= 150),
            midpoint   = (frag_start + frag_end) // 2,
            mapq       = read.mapping_quality,
        )