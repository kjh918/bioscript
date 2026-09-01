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
    [MODIFIED] 빈(Bin)을 삭제(Drop)하지 않고, 필터링 사유를 Column으로 태깅합니다.
    유전체의 물리적 연속성을 유지하여 추후 Segmentation 알고리즘의 정확도를 보장합니다.
    """
    log(f"Applying final filters (Tagging)... (Input bins: {len(df)})")

    df_final = df.copy()

    # 1. 조건별 Boolean Mask 생성
    # GC 함량이 정상 범위를 벗어났는가?
    mask_abnormal_gc = ~df_final['gc'].between(args.MinGC, args.MaxGC)
    
    # Mappability가 기준치 미달인가?
    mask_low_map = df_final['mappability'] < args.MinMappability
    
    # Blacklist 포함 여부 (기존 컬럼이 없다면 False로 초기화)
    if 'is_blacklisted' in df_final.columns:
        mask_blacklist = df_final['is_blacklisted'] == True
    else:
        mask_blacklist = pd.Series(False, index=df_final.index)

    # 2. 내부 연산용 Boolean Column 추가 (True/False가 Pandas 연산 속도에 압도적으로 유리함)
    df_final['is_abnormal_gc'] = mask_abnormal_gc
    df_final['is_low_mappability'] = mask_low_map
    df_final['is_blacklisted'] = mask_blacklist
    
    # 세 가지 중 하나라도 걸리면 최종적으로 필터링(사용 불가) 됨을 의미
    df_final['is_filtered'] = mask_abnormal_gc | mask_low_map | mask_blacklist

    # 3. [선택] 연구원님 가독성/리포트 출력을 위한 'O' / 'X' 문자열 매핑 컬럼 추가
    # 연산 시에는 Boolean을 쓰고, 엑셀이나 눈으로 볼 때는 이 컬럼을 보시면 됩니다.
    df_final['low_map_flag'] = np.where(df_final['is_low_mappability'], 'O', 'X')
    df_final['filtered_flag'] = np.where(df_final['is_filtered'], 'O', 'X')

    dropped_count = df_final['is_filtered'].sum()
    log(f"Tagging complete. (Total bins: {len(df_final)}, Flagged as filtered: {dropped_count})")
    
    # 행을 삭제하지 않고 원본 길이 그대로 반환!
    return df_final