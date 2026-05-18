import pandas as pd
import numpy as np
import pysam
import pyBigWig
from pyfaidx import Fasta
from utils import log

def get_chromosomes(fasta_path, include_sex=True):
    fa = pysam.FastaFile(fasta_path)
    chroms = [c for c in fa.references if c not in ["chrM", "MT", "M"]]
    if not include_sex:
        chroms = [c for c in chroms if c not in ["chrX", "chrY", "X", "Y"]]
    fa.close()
    return chroms

def generate_bins(fasta_path, bin_size, chromosomes):
    fa = pysam.FastaFile(fasta_path)
    rows = []
    for chrom in chromosomes:
        c_len = fa.get_reference_length(chrom)
        for start in range(0, c_len, bin_size):
            end = min(start + bin_size, c_len)
            rows.append((chrom, start, end))
    fa.close()
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])

def annotate_bin_metadata(df, fasta_path, mappability_bw=None):
    """
    [통합 함수] GC, 유효 길이, Mappability를 한 번의 루프로 주석 처리합니다.
    """
    fa = Fasta(fasta_path)
    bw = pyBigWig.open(mappability_bw) if mappability_bw else None
    
    results = []
    for r in df.itertuples():
        seq = fa[r.chrom][r.start:r.end].seq.upper()
        eff_len = sum(base in "ACGT" for base in seq)
        gc = (seq.count("G") + seq.count("C")) / eff_len if eff_len > 0 else np.nan
        mappability = bw.stats(r.chrom, r.start, r.end, type="mean")[0] if bw else 1.0
        mappability = mappability if mappability is not None else 0.0
        mappability = float(np.clip(mappability, 0.0, 1.0))
        results.append({
            'gc': gc,
            'mappability': mappability
        })
    
    if bw: bw.close()
    return pd.concat([df, pd.DataFrame(results)], axis=1)


def apply_final_filters(df, args):
    """
    모든 정적/동적 지표를 결합하여 분석에 사용할 최종 빈(Bin)을 선별합니다.
    """
    log(f"Applying final filters... (Input bins: {len(df)})")

    # 1. [정적 필터] 유전체 구조 기준
    # GC 함량 범위 및 Mappability 기준 통과 여부
    mask = (df['gc'].between(args.MinGC, args.MaxGC))
    mask &= (df['mappability'] >= args.MinMappability)
    
    # Blacklist BED 파일이 적용되었다면 해당 구간 제외
    if 'is_blacklisted' in df.columns:
        mask &= (df['is_blacklisted'] == False)

    df_final = df[mask].copy().reset_index(drop=True)
    
    log(f"Filtering complete. (Output bins: {len(df_final)}, Dropped: {len(df) - len(df_final)})")
    
    return df_final