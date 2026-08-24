"""
wps/peak_calling.py
===================
Snyder et al. 2016 방법론 기반 nucleosome peak calling

입력: L-WPS norm.npy (median-subtracted + SG-smoothed)
출력: peaks BED (chrom, start, end, center, score, region_len)

알고리즘
--------
1. above-zero segmentation
   - 연속 양수 구간 추출 (최대 5bp gap 허용)

2. region 크기별 처리
   - <50bp  또는 >450bp : discard
   - 50~150bp  : median 위 max-sum contiguous window → peak 1개
   - 150~450bp : 위 동일 로직으로 50~150bp sub-window 모두 report

3. score
   - max(window) - mean(left_min, right_min)
   - 150~450bp multi-window: neighboring minima = 0

출력 컬럼
---------
chrom, start, end, center, score, region_len, window_len, source_region
"""
from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .utils import HG38_CHROM_SIZES, STANDARD_CHROMS, read_manifest, load_track

log = logging.getLogger(__name__)

# ── 상수 ──────────────────────────────────────────────────────────────
_GAP_ALLOW   = 5    # above-zero segmentation 허용 gap (bp)
_REG_MIN     = 50   # region 최소 크기
_REG_MID     = 150  # single / multi 경계
_REG_MAX     = 450  # region 최대 크기
_WIN_MIN     = 50   # sub-window 최소 크기
_WIN_MAX     = 150  # sub-window 최대 크기


# ── 데이터 클래스 ──────────────────────────────────────────────────────
@dataclass
class Peak:
    chrom:         str
    start:         int
    end:           int
    center:        int
    score:         float
    region_len:    int
    window_len:    int
    source_region: str   # "50-150" | "150-450"


# ── above-zero segmentation ───────────────────────────────────────────
def _segment_above_zero(
    arr: np.ndarray,
    gap_allow: int = _GAP_ALLOW,
) -> list[tuple[int, int]]:
    """
    arr에서 0 초과 구간을 추출.
    gap_allow 이하의 연속 0-이하 구간은 무시(bridge).

    Returns
    -------
    [(start, end), ...]  0-based half-open 좌표
    """
    above = (arr > 0).astype(np.int8)
    n     = len(above)
    segs: list[tuple[int, int]] = []

    i = 0
    while i < n:
        if above[i] == 0:
            i += 1
            continue

        # region 시작
        seg_start = i
        j = i + 1
        gap = 0

        while j < n:
            if above[j] > 0:
                gap = 0
            else:
                gap += 1
                if gap > gap_allow:
                    break
            j += 1

        seg_end = j - gap  # gap 앞까지
        segs.append((seg_start, seg_end))
        i = j

    return segs


# ── max-sum contiguous window ─────────────────────────────────────────
def _max_sum_window(
    arr:      np.ndarray,
    win_min:  int,
    win_max:  int,
    median_threshold: float = 0.0,
) -> Optional[tuple[int, int, float]]:
    """
    arr[median_threshold 초과 위치]에서 길이 win_min~win_max 인
    최대 합 연속 window를 찾아 (start, end, window_sum) 반환.
    못 찾으면 None.
    """
    above_med = np.where(arr > median_threshold, arr, 0.0)
    n         = len(above_med)

    best_sum   = -np.inf
    best_start = -1
    best_end   = -1

    for w in range(win_min, min(win_max + 1, n + 1)):
        # sliding window sum
        cumsum = np.concatenate([[0], np.cumsum(above_med)])
        sums   = cumsum[w:] - cumsum[:n - w + 1]
        if sums.size == 0:
            continue
        idx = int(np.argmax(sums))
        if sums[idx] > best_sum:
            best_sum   = float(sums[idx])
            best_start = idx
            best_end   = idx + w

    if best_start < 0:
        return None
    return best_start, best_end, best_sum


# ── 인접 minima 찾기 ──────────────────────────────────────────────────
def _adjacent_minima(
    full_arr:    np.ndarray,
    reg_start:   int,
    reg_end:     int,
    search_bp:   int = 300,
) -> tuple[float, float]:
    """
    region 왼쪽/오른쪽 search_bp 범위에서 최솟값 반환.
    범위 밖이면 0.0.
    """
    n = len(full_arr)

    left_start = max(0, reg_start - search_bp)
    left_end   = reg_start
    right_start = reg_end
    right_end   = min(n, reg_end + search_bp)

    left_min  = float(np.min(full_arr[left_start:left_end]))  if left_end  > left_start  else 0.0
    right_min = float(np.min(full_arr[right_start:right_end])) if right_end > right_start else 0.0

    return left_min, right_min


# ── 단일 region → peak 리스트 ─────────────────────────────────────────
def _call_peaks_from_region(
    chrom:     str,
    arr:       np.ndarray,   # region 전체 (full chrom 배열)
    reg_start: int,
    reg_end:   int,
    offset:    int = 0,      # chrom 내 절대 좌표 변환용
) -> list[Peak]:
    """
    단일 above-zero region에서 peak 추출.
    """
    reg_len = reg_end - reg_start
    seg     = arr[reg_start:reg_end]

    if reg_len < _REG_MIN or reg_len > _REG_MAX:
        return []

    peaks: list[Peak] = []

    if _REG_MIN <= reg_len <= _REG_MID:
        # ── 50~150bp: single peak ───────────────────────────────────
        med    = float(np.median(seg))
        result = _max_sum_window(seg, _WIN_MIN, _WIN_MAX, median_threshold=med)
        if result is None:
            return []

        w_start, w_end, _ = result
        w_arr = seg[w_start:w_end]
        score_max = float(np.max(w_arr))

        left_min, right_min = _adjacent_minima(arr, reg_start, reg_end)
        score = score_max - (left_min + right_min) / 2.0

        abs_start  = offset + reg_start + w_start
        abs_end    = offset + reg_start + w_end
        abs_center = (abs_start + abs_end) // 2

        peaks.append(Peak(
            chrom         = chrom,
            start         = abs_start,
            end           = abs_end,
            center        = abs_center,
            score         = round(score, 4),
            region_len    = reg_len,
            window_len    = w_end - w_start,
            source_region = "50-150",
        ))

    elif _REG_MID < reg_len <= _REG_MAX:
        # ── 150~450bp: 50~150bp sub-window 모두 report ──────────────
        med = float(np.median(seg))

        # sliding으로 모든 50~150bp window를 탐색
        for w in range(_WIN_MIN, _WIN_MAX + 1):
            if w > len(seg):
                break
            cumsum = np.concatenate([[0], np.cumsum(
                np.where(seg > med, seg, 0.0)
            )])
            sums = cumsum[w:] - cumsum[:len(seg) - w + 1]
            if sums.size == 0:
                continue

            # 겹치지 않는 최대 window를 탐욕적으로 선택
            used = np.zeros(len(seg), dtype=bool)
            order = np.argsort(-sums)

            for idx in order:
                if sums[idx] <= 0:
                    break
                if used[idx:idx + w].any():
                    continue
                used[idx:idx + w] = True

                w_arr     = seg[idx:idx + w]
                score_max = float(np.max(w_arr))
                # 150~450 bp: neighboring minima = 0
                score = score_max - 0.0

                abs_start  = offset + reg_start + idx
                abs_end    = offset + reg_start + idx + w
                abs_center = (abs_start + abs_end) // 2

                peaks.append(Peak(
                    chrom         = chrom,
                    start         = abs_start,
                    end           = abs_end,
                    center        = abs_center,
                    score         = round(score, 4),
                    region_len    = reg_len,
                    window_len    = w,
                    source_region = "150-450",
                ))
            break  # w 루프: 대표 window 크기 1개만 (논문 방식)

    return peaks


# ── 염색체 단위 peak calling ──────────────────────────────────────────
def call_peaks_chrom(
    chrom:    str,
    norm_arr: np.ndarray,
    offset:   int = 0,
) -> list[Peak]:
    """
    단일 염색체 norm array → Peak 리스트.

    Parameters
    ----------
    chrom    : 염색체 이름
    norm_arr : L-WPS normalized (median-subtracted + SG-smoothed)
    offset   : 절대 좌표 변환 (bin 단위 slicing 시 사용)
    """
    segs  = _segment_above_zero(norm_arr, _GAP_ALLOW)
    peaks: list[Peak] = []

    for reg_start, reg_end in segs:
        peaks.extend(
            _call_peaks_from_region(chrom, norm_arr, reg_start, reg_end, offset)
        )

    return peaks


# ── BED 저장 ──────────────────────────────────────────────────────────
def write_peaks_bed(peaks: list[Peak], output_path: str) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tcenter\tscore\t"
            "region_len\twindow_len\tsource_region\n"
        )
        for p in peaks:
            fh.write(
                f"{p.chrom}\t{p.start}\t{p.end}\t{p.center}\t"
                f"{p.score}\t{p.region_len}\t{p.window_len}\t"
                f"{p.source_region}\n"
            )


# ── multiprocessing worker ────────────────────────────────────────────
def _call_peaks_chrom_worker(
    chrom: str,
    wps_dir: str,
    prefix: str,
    mode: str,
) -> tuple[str, list[Peak], str | None]:
    """
    ProcessPoolExecutor에서 실행되는 chromosome 단위 worker.

    각 process가 자신의 norm.npy를 직접 로드하므로,
    큰 chromosome array를 parent process에서 pickle하여 전달하지 않는다.
    """
    try:
        manifest = read_manifest(wps_dir, mode)

        arr = load_track(
            wps_dir=wps_dir,
            mode=mode,
            chrom=chrom,
            metric="norm",
            prefix=prefix,
            manifest=manifest,
        )

        if arr is None:
            return chrom, [], f"norm.npy 없음: {chrom} — skip"

        # 가능한 경우 불필요한 copy 방지
        arr = np.asarray(arr, dtype=np.float32)
        arr = np.nan_to_num(arr, nan=0.0)

        peaks = call_peaks_chrom(chrom, arr)
        return chrom, peaks, None

    except Exception as exc:
        return chrom, [], f"{type(exc).__name__}: {exc}"


# ── 공개 인터페이스 ────────────────────────────────────────────────────
def run(
    wps_dir:    str,
    output_bed: str,
    prefix:     str = "",
    chroms:     Optional[list[str]] = None,
    mode:       str = "L",
    n_jobs:     Optional[int] = None,
) -> list[Peak]:
    """
    manifest NPY에서 L-WPS norm 로드 → chromosome 병렬 peak calling → BED 저장

    Parameters
    ----------
    wps_dir    : compute.run() out_dir 루트
    output_bed : 출력 BED 경로
    prefix     : sample_id prefix
    chroms     : 처리할 염색체 (None이면 standard 전체)
    mode       : "L" 고정 (S-WPS peak calling은 미지원)
    n_jobs     : 동시에 처리할 chromosome process 수.
                 None이면 min(CPU 수, chromosome 수) 사용.
                 메모리가 부족하면 2~4 정도로 제한 권장.

    Returns
    -------
    전체 Peak 리스트
    """
    if mode != "L":
        raise ValueError("peak calling은 L-WPS 전용입니다.")

    target = list(chroms or STANDARD_CHROMS)
    if not target:
        write_peaks_bed([], output_bed)
        return []

    if n_jobs is None:
        n_jobs = min(len(target), os.cpu_count() or 1)
    else:
        n_jobs = max(1, min(int(n_jobs), len(target)))

    log.info(
        "peak calling 시작: chromosomes=%d, n_jobs=%d",
        len(target),
        n_jobs,
    )

    # 결과를 chromosome별로 저장한 뒤 target 순서대로 합친다.
    # as_completed 순서로 바로 append하면 BED chromosome 순서가 뒤섞일 수 있다.
    peaks_by_chrom: dict[str, list[Peak]] = {chrom: [] for chrom in target}

    if n_jobs == 1:
        # 디버깅/저메모리 환경용 serial fallback
        for chrom in target:
            chrom_name, chrom_peaks, error = _call_peaks_chrom_worker(
                chrom=chrom,
                wps_dir=wps_dir,
                prefix=prefix,
                mode=mode,
            )

            if error:
                log.warning("%s", error)
                continue

            peaks_by_chrom[chrom_name] = chrom_peaks
            log.info("  %s: %d peaks", chrom_name, len(chrom_peaks))

    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            future_to_chrom = {
                executor.submit(
                    _call_peaks_chrom_worker,
                    chrom,
                    wps_dir,
                    prefix,
                    mode,
                ): chrom
                for chrom in target
            }

            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]

                try:
                    chrom_name, chrom_peaks, error = future.result()
                except Exception:
                    log.exception("  %s: worker 실패", chrom)
                    continue

                if error:
                    log.warning("%s", error)
                    continue

                peaks_by_chrom[chrom_name] = chrom_peaks
                log.info("  %s: %d peaks", chrom_name, len(chrom_peaks))

    all_peaks: list[Peak] = []
    for chrom in target:
        all_peaks.extend(peaks_by_chrom.get(chrom, []))

    write_peaks_bed(all_peaks, output_bed)
    log.info(
        "peak calling 완료: %d peaks → %s",
        len(all_peaks),
        output_bed,
    )
    return all_peaks
