import pandas as pd
import numpy as np

def calculate_titv_ratio(df, ref_col="ref", alt_col="alt"):
    """
    DataFrame에서 Reference와 Alternate 염기 정보를 받아 
    Transition / Transversion (Ti/Tv) 비율을 계산합니다.
    
    Parameters:
    -----------
    df : pd.DataFrame
        돌연변이 정보가 담긴 데이터프레임
    ref_col : str
        Reference 염기가 있는 컬럼명 (기본값: 'ref')
    alt_col : str
        Alternate 염기가 있는 컬럼명 (기본값: 'alt')
        
    Returns:
    --------
    dict
        ti_count (Transition 개수), 
        tv_count (Transversion 개수), 
        titv_ratio (Ti/Tv 비율)을 포함한 딕셔너리
    """
    if df is None or df.empty or ref_col not in df.columns or alt_col not in df.columns:
        return {"ti_count": 0, "tv_count": 0, "titv_ratio": 0.0}

    # 1. 단일 염기 다형성(SNP)만 필터링 (Indel 제외)
    valid_bases = ['A', 'C', 'G', 'T']
    
    # 대문자로 변환
    refs = df[ref_col].astype(str).str.upper()
    alts = df[alt_col].astype(str).str.upper()
    
    # SNP 마스크 생성 (ref와 alt가 모두 A, C, G, T 중 하나이고 서로 다른 경우)
    snp_mask = refs.isin(valid_bases) & alts.isin(valid_bases) & (refs != alts)
    
    filtered_refs = refs[snp_mask]
    filtered_alts = alts[snp_mask]
    
    # 2. Transition (Ti) 정의: 퓨린(A, G) <-> 퓨린, 피리미딘(C, T) <-> 피리미딘
    # A <-> G, C <-> T
    transitions = {
        ('A', 'G'), ('G', 'A'),
        ('C', 'T'), ('T', 'C')
    }
    
    # 3. Ti / Tv 카운팅
    ti_count = 0
    tv_count = 0
    
    for r, a in zip(filtered_refs, filtered_alts):
        if (r, a) in transitions:
            ti_count += 1
        else:
            tv_count += 1
            
    # 4. 비율 계산 (ZeroDivision 방지)
    titv_ratio = float(ti_count / tv_count) if tv_count > 0 else float('inf')
    
    return {
        "ti_count": ti_count,
        "tv_count": tv_count,
        "titv_ratio": titv_ratio
    }

# =====================================================================
# 사용 예시 (테스트 코드)
# =====================================================================
if __name__ == "__main__":
    # 가상의 VCF 파싱 데이터 생성
    data = [
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}, # Ti
        {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T"}, # Ti
        {"chrom": "chr1", "pos": 300, "ref": "A", "alt": "C"}, # Tv
        {"chrom": "chr1", "pos": 400, "ref": "G", "alt": "T"}, # Tv
        {"chrom": "chr1", "pos": 500, "ref": "T", "alt": "C"}, # Ti
        {"chrom": "chr1", "pos": 600, "ref": "A", "alt": "AT"},# Indel (무시됨)
        {"chrom": "chr1", "pos": 700, "ref": "G", "alt": "A"}, # Ti
    ]
    
    df_variants = pd.DataFrame(data)
    
    result = calculate_titv_ratio(df_variants, ref_col="ref", alt_col="alt")
    
    print("=== Transition / Transversion (Ti/Tv) Analysis ===")
    print(f"Total valid SNPs: {result['ti_count'] + result['tv_count']}")
    print(f"Transitions (Ti) : {result['ti_count']}")
    print(f"Transversions (Tv): {result['tv_count']}")
    print(f"Ti/Tv Ratio      : {result['titv_ratio']:.3f}")