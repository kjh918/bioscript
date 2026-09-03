"""
nipt_loader.py
==============
실제 cbNIPT 파이프라인 output TSV 로더.

Input 1: cnv.tsv (전체 genome bin)
  columns: chrom, start, end, bin_id, bin_BAF, copy_number_signal

Input 2: nipt_results.tsv (syndrome 판정 결과)
  columns: SYNDROME, NIPT_GROUP, FEAT_RANK, FEATURE_TYPE, FEATURE_NAME,
           CHROMOSOME, OVERLAP_BINS, OBSERVED_CN, BAF_MEDIAN,
           HETERO_RATE, HOMO_RATE, DETECTED_CNV, DIAGNOSIS, EVIDENCE,
           LOW_RESOLUTION_WARNING
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
import pandas as pd

from .nipt_markers import NiptSyndrome, MarkerFeature, CALL_COLORS
from .models import SampleData, CnvEvent

# ---------------------------------------------------------------------------
# HG38 Chromosome position Map
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# GRCh38 chromosome map
# - chromosome length: GRCh38 primary assembly
# - centromere: NCBI GRC modeled centromere
# - centromere boundaries expanded to 100 kb bin boundaries
#
# 좌표계:
#   0-based half-open [start, end)
#
# 예:
#   chr1 first bin = [0, 100000)
# ---------------------------------------------------------------------------

_CHR_MAP = {
    "assembly": "GRCh38",
    "version": "GRCh38.p14",
    "bin_size": 100_000,
    "source": "NCBI Genome Reference Consortium GRCh38.p14",
    "chromosomes": {
        "1": {
            "length": 248_956_422,
            "p_arm":      {"start": 0,           "end": 122_000_000},
            "centromere": {"start": 122_000_000, "end": 125_200_000},
            "q_arm":      {"start": 125_200_000, "end": 248_956_422},
        },
        "2": {
            "length": 242_193_529,
            "p_arm":      {"start": 0,          "end": 92_100_000},
            "centromere": {"start": 92_100_000, "end": 94_100_000},
            "q_arm":      {"start": 94_100_000, "end": 242_193_529},
        },
        "3": {
            "length": 198_295_559,
            "p_arm":      {"start": 0,          "end": 90_700_000},
            "centromere": {"start": 90_700_000, "end": 93_700_000},
            "q_arm":      {"start": 93_700_000, "end": 198_295_559},
        },
        "4": {
            "length": 190_214_555,
            "p_arm":      {"start": 0,          "end": 49_700_000},
            "centromere": {"start": 49_700_000, "end": 51_800_000},
            "q_arm":      {"start": 51_800_000, "end": 190_214_555},
        },
        "5": {
            "length": 181_538_259,
            "p_arm":      {"start": 0,          "end": 46_400_000},
            "centromere": {"start": 46_400_000, "end": 50_100_000},
            "q_arm":      {"start": 50_100_000, "end": 181_538_259},
        },
        "6": {
            "length": 170_805_979,
            "p_arm":      {"start": 0,          "end": 58_500_000},
            "centromere": {"start": 58_500_000, "end": 59_900_000},
            "q_arm":      {"start": 59_900_000, "end": 170_805_979},
        },
        "7": {
            "length": 159_345_973,
            "p_arm":      {"start": 0,          "end": 58_100_000},
            "centromere": {"start": 58_100_000, "end": 60_900_000},
            "q_arm":      {"start": 60_900_000, "end": 159_345_973},
        },
        "8": {
            "length": 145_138_636,
            "p_arm":      {"start": 0,          "end": 44_000_000},
            "centromere": {"start": 44_000_000, "end": 45_900_000},
            "q_arm":      {"start": 45_900_000, "end": 145_138_636},
        },
        "9": {
            "length": 138_394_717,
            "p_arm":      {"start": 0,          "end": 43_200_000},
            "centromere": {"start": 43_200_000, "end": 45_600_000},
            "q_arm":      {"start": 45_600_000, "end": 138_394_717},
        },
        "10": {
            "length": 133_797_422,
            "p_arm":      {"start": 0,          "end": 39_600_000},
            "centromere": {"start": 39_600_000, "end": 41_600_000},
            "q_arm":      {"start": 41_600_000, "end": 133_797_422},
        },
        "11": {
            "length": 135_086_622,
            "p_arm":      {"start": 0,          "end": 51_000_000},
            "centromere": {"start": 51_000_000, "end": 54_500_000},
            "q_arm":      {"start": 54_500_000, "end": 135_086_622},
        },
        "12": {
            "length": 133_275_309,
            "p_arm":      {"start": 0,          "end": 34_700_000},
            "centromere": {"start": 34_700_000, "end": 37_200_000},
            "q_arm":      {"start": 37_200_000, "end": 133_275_309},
        },
        "13": {
            "length": 114_364_328,
            "p_arm":      {"start": 0,          "end": 16_000_000},
            "centromere": {"start": 16_000_000, "end": 18_100_000},
            "q_arm":      {"start": 18_100_000, "end": 114_364_328},
        },
        "14": {
            "length": 107_043_718,
            "p_arm":      {"start": 0,          "end": 16_000_000},
            "centromere": {"start": 16_000_000, "end": 18_200_000},
            "q_arm":      {"start": 18_200_000, "end": 107_043_718},
        },
        "15": {
            "length": 101_991_189,
            "p_arm":      {"start": 0,          "end": 17_000_000},
            "centromere": {"start": 17_000_000, "end": 19_800_000},
            "q_arm":      {"start": 19_800_000, "end": 101_991_189},
        },
        "16": {
            "length": 90_338_345,
            "p_arm":      {"start": 0,          "end": 36_300_000},
            "centromere": {"start": 36_300_000, "end": 38_300_000},
            "q_arm":      {"start": 38_300_000, "end": 90_338_345},
        },
        "17": {
            "length": 83_257_441,
            "p_arm":      {"start": 0,          "end": 22_800_000},
            "centromere": {"start": 22_800_000, "end": 26_900_000},
            "q_arm":      {"start": 26_900_000, "end": 83_257_441},
        },
        "18": {
            "length": 80_373_285,
            "p_arm":      {"start": 0,          "end": 15_400_000},
            "centromere": {"start": 15_400_000, "end": 20_900_000},
            "q_arm":      {"start": 20_900_000, "end": 80_373_285},
        },
        "19": {
            "length": 58_617_616,
            "p_arm":      {"start": 0,          "end": 24_400_000},
            "centromere": {"start": 24_400_000, "end": 27_200_000},
            "q_arm":      {"start": 27_200_000, "end": 58_617_616},
        },
        "20": {
            "length": 64_444_167,
            "p_arm":      {"start": 0,          "end": 26_400_000},
            "centromere": {"start": 26_400_000, "end": 30_100_000},
            "q_arm":      {"start": 30_100_000, "end": 64_444_167},
        },
        "21": {
            "length": 46_709_983,
            "p_arm":      {"start": 0,          "end": 10_800_000},
            "centromere": {"start": 10_800_000, "end": 13_000_000},
            "q_arm":      {"start": 13_000_000, "end": 46_709_983},
        },
        "22": {
            "length": 50_818_468,
            "p_arm":      {"start": 0,          "end": 12_900_000},
            "centromere": {"start": 12_900_000, "end": 15_100_000},
            "q_arm":      {"start": 15_100_000, "end": 50_818_468},
        },
        "X": {
            "length": 156_040_895,
            "p_arm":      {"start": 0,          "end": 58_600_000},
            "centromere": {"start": 58_600_000, "end": 62_500_000},
            "q_arm":      {"start": 62_500_000, "end": 156_040_895},
        },
        "Y": {
            "length": 57_227_415,
            "p_arm":      {"start": 0,          "end": 10_300_000},
            "centromere": {"start": 10_300_000, "end": 10_600_000},
            "q_arm":      {"start": 10_600_000, "end": 57_227_415},
        },
    },
}
# ---------------------------------------------------------------------------
# DIAGNOSIS → call 변환
# ---------------------------------------------------------------------------
_DIAG_MAP = {
    # 확실한 이상
    "GAIN":             "HIGH_RISK",
    "LOSS":             "HIGH_RISK",
    "AMPLIFICATION":    "HIGH_RISK",
    "DELETION":         "HIGH_RISK",
    "MOSAIC_GAIN":      "SUSPECTED",
    "MOSAIC_LOSS":      "SUSPECTED",
    "SUSPECTED_GAIN":   "SUSPECTED",
    "SUSPECTED_LOSS":   "SUSPECTED",
    "SUSPECTED":        "SUSPECTED",
    # 정상
    "NEUT":             "LOW_RISK",
    "NEGATIVE":         "LOW_RISK",
    "NORMAL":           "LOW_RISK",
    "NO_CALL":          "LOW_RISK",
}

_DETECTED_CNV_MAP = {
    "GAIN":   "HIGH_RISK",
    "LOSS":   "HIGH_RISK",
    "NEUT":   "LOW_RISK",
    "AMP":    "HIGH_RISK",
    "DEL":    "HIGH_RISK",
}


def _to_call(diagnosis: str, detected_cnv: str = "") -> str:
    d = str(diagnosis).strip().upper()
    c = str(detected_cnv).strip().upper()
    return (
        _DIAG_MAP.get(d)
        or _DETECTED_CNV_MAP.get(c)
        or "LOW_RISK"
    )


# ---------------------------------------------------------------------------
# CNV TSV loader (전체 genome)
# ---------------------------------------------------------------------------
_CNV_COL_MAP = {
    # 실제 파이프라인 컬럼 → 표준 컬럼
    "copy_number_signal": "cn",
    "bin_baf":            "baf",
    "bin_BAF":            "baf",
    "start":              "pos",
}


def load_cnv_tsv(path: "str | Path") -> dict[str, pd.DataFrame]:
    """
    전체 genome cnv.tsv → {chrom: DataFrame(pos, cn, baf, start, end)}
    chrom 컬럼 기준으로 분리.
    """
    df = pd.read_csv(Path(path), sep="\t")
    df.columns = [c.strip() for c in df.columns]

    # 컬럼 정규화
    rename = {}
    for col in df.columns:
        low = col.lower()
        if low in _CNV_COL_MAP:
            rename[col] = _CNV_COL_MAP[low]
        elif col in _CNV_COL_MAP:
            rename[col] = _CNV_COL_MAP[col]
    df = df.rename(columns=rename)

    # pos 없으면 start 사용
    if "pos" not in df.columns and "start" in df.columns:
        df["pos"] = df["start"]

    # chrom 정규화
    df["chrom"] = df["chrom"].str.replace("chr", "", regex=False).str.upper()

    result = {}
    for chrom, grp in df.groupby("chrom"):
        sub = grp.reset_index(drop=True)
        result[str(chrom)] = sub

    print(f"[nipt_loader] cnv: {len(df):,} bins, {len(result)} chroms")
    return result


# ---------------------------------------------------------------------------
# NIPT results TSV loader
# ---------------------------------------------------------------------------
def load_nipt_results(
    path: "str | Path",
    marker_path: "str | Path | None" = None,
) -> dict[str, NiptSyndrome]:
    """
    nipt_results.tsv → {syndrome_key: NiptSyndrome}
    marker_path: nipt_markers.tsv — feature 좌표 조회용
    """
    # ── marker_path에서 feature 좌표 미리 로드 ────────────────────────────
    # key: (syndrome_upper, feat_name) → (start, end)
    feat_coords: dict[tuple, tuple] = {}
    if marker_path and Path(marker_path).exists():
        mdf = pd.read_csv(Path(marker_path), sep="\t")
        mdf.columns = [c.strip() for c in mdf.columns]
        for _, mrow in mdf.iterrows():
            syn_key  = str(mrow["SYNDROME"]).strip().upper().replace(" ", "_").replace("-", "_")
            feat_key = str(mrow["FEATURE_NAME"]).strip()
            start    = int(mrow.get("GENOMIC_POS_START", 1))
            end      = int(mrow.get("GENOMIC_POS_END", 250_000_000))
            feat_coords[(syn_key, feat_key)] = (start, end)
            # nipt_id 기반 키도 등록
            nipt_id_key = str(mrow.get("NIPT_ID", "")).strip()
            if nipt_id_key:
                feat_coords[(nipt_id_key, feat_key)] = (start, end)
        print(f"[nipt_loader] marker coords loaded: {len(feat_coords)} entries")

    df = pd.read_csv(Path(path), sep="\t")
    df.columns = [c.strip() for c in df.columns]

    col_map = {c.upper(): c for c in df.columns}
    def _col(name: str):
        return col_map.get(name.upper(), name)

    syndromes: dict[str, NiptSyndrome] = {}
    _rank = {"HIGH_RISK": 3, "SUSPECTED": 2, "LOW_RISK": 1, "UNKNOWN": 0}

    for _, row in df.iterrows():
        syndrome  = str(row[_col("SYNDROME")]).strip()
        group     = str(row[_col("NIPT_GROUP")]).strip()
        feat_rank = int(row.get(_col("FEAT_RANK"), 9))
        feat_type = str(row[_col("FEATURE_TYPE")]).strip()
        feat_name = str(row[_col("FEATURE_NAME")]).strip()
        chrom     = str(row[_col("CHROMOSOME")]).strip()
        diagnosis = str(row.get(_col("DIAGNOSIS"), "")).strip()
        detected  = str(row.get(_col("DETECTED_CNV"), "")).strip()
        observed_cn = row.get(_col("OBSERVED_CN"))
        observed_cn = float(observed_cn) if pd.notna(observed_cn) else None
        baf_median  = row.get(_col("BAF_MEDIAN"))
        baf_median  = float(baf_median) if pd.notna(baf_median) else None

        call    = _to_call(diagnosis, detected)
        nipt_id = syndrome.upper().replace(" ", "_").replace("-", "_")

        if nipt_id not in syndromes:
            syndromes[nipt_id] = NiptSyndrome(
                nipt_id  = nipt_id,
                syndrome = syndrome,
                group    = group,
                features = [],
                call     = "LOW_RISK",
                cn_value = None,
            )

        syn = syndromes[nipt_id]

        if _rank.get(call, 0) > _rank.get(syn.call, 0):
            syn.call = call

        if feat_rank <= 3 and observed_cn is not None:
            syn.cn_value   = observed_cn
            syn.baf_median = baf_median

        # feature 좌표 — marker_path 우선, 없으면 chrom 전체
        existing = {f.feature_name for f in syn.features}
        if feat_name not in existing:
            coords = (
                feat_coords.get((nipt_id, feat_name))
                or feat_coords.get((nipt_id, chrom))
                or feat_coords.get((syndrome.upper().replace(" ", "_"), feat_name))
            )
            if coords:
                start, end = coords
            else:
                # fallback: chrom 전체
                from .reference import CHROM_SIZES
                chrom_key = chrom.replace("chr", "").upper()
                end   = CHROM_SIZES.get(chrom_key, 250_000_000)
                start = 1

            syn.features.append(MarkerFeature(
                feature_name = feat_name,
                feature_type = feat_type,
                chromosome   = chrom,
                start        = start,
                end          = end,
            ))

    print(f"[nipt_loader] results: {len(syndromes)} syndromes")
    _log_calls(syndromes)
    return syndromes


def _log_calls(syndromes: dict[str, NiptSyndrome]):
    for s in syndromes.values():
        if s.call != "LOW_RISK":
            cn = f"CN={s.cn_value:.3f}" if s.cn_value is not None else ""
            print(f"  [{s.call}] {s.syndrome} {cn}")


# ---------------------------------------------------------------------------
# CnvEvent 자동 생성 (HIGH_RISK / SUSPECTED syndrome → events)
# ---------------------------------------------------------------------------
_TYPE_MAP = {
    "Autosome Abnormality": {
        "HIGH_RISK": lambda cn: "trisomy"  if cn > 2.5 else "monosomy",
        "SUSPECTED": lambda cn: "suspected_gain" if cn > 2 else "suspected_loss",
    },
    "Sex Chromosome Abnormality": {
        "HIGH_RISK": lambda cn: "gain" if cn > 2.5 else "monosomy",
        "SUSPECTED": lambda cn: "suspected_gain" if cn > 2 else "suspected_loss",
    },
    "Micro Deletion": {
        "HIGH_RISK": lambda cn: "partial_loss",
        "SUSPECTED": lambda cn: "suspected_loss",
    },
}

_CALL_COLORS = {
    "HIGH_RISK": "#E53E3E",
    "SUSPECTED": "#DD6B20",
    "LOW_RISK":  "#38A169",
}

_ISCN_PREFIX = {
    "trisomy":        "+",
    "monosomy":       "-",
    "gain":           "dup(",
    "partial_loss":   "del(",
    "suspected_gain": "dup(",
    "suspected_loss": "del(",
}


def build_events_from_results(
    syndromes: dict[str, NiptSyndrome],
) -> list[CnvEvent]:
    """HIGH_RISK / SUSPECTED syndrome에서 CnvEvent 자동 생성."""
    events = []
    seen = set()

    for syn in syndromes.values():
        if syn.call not in ("HIGH_RISK", "SUSPECTED"):
            continue

        primary = next(
            (f for f in syn.features
             if f.feature_type in ("TargetChromosome", "PrimaryTargetRegion")),
            syn.features[0] if syn.features else None,
        )
        if primary is None:
            continue

        chrom = primary.chromosome.replace("chr", "")
        key   = (chrom, syn.call)
        if key in seen:
            continue
        seen.add(key)

        cn = syn.cn_value or 2.0
        type_fn = _TYPE_MAP.get(syn.group, {}).get(syn.call)
        ev_type = type_fn(cn) if type_fn else "unknown"

        # ISCN 표기
        prefix = _ISCN_PREFIX.get(ev_type, "?")
        if ev_type in ("trisomy", "monosomy"):
            iscn_str = prefix + chrom
        else:
            feat = primary.feature_name
            iscn_str = f"{prefix}{feat})"

        events.append(CnvEvent(
            chr   = chrom,
            type  = ev_type,
            cn    = round(cn, 2),
            iscn  = iscn_str,
            start = primary.start,
            stop  = primary.end,
            color = _CALL_COLORS.get(syn.call, "#718096"),
        ))

    return events
