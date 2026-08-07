from pathlib import Path
import pandas as pd


def _read_table(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix in {'.tsv', '.txt', '.bed'}:
        return pd.read_csv(path, sep='\t')
    if suffix == '.csv':
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine='python')


def _norm_chrom(series):
    return series.astype(str).str.replace('chr', '', regex=False).str.upper()


def load_signal(path):
    """
    Required: chrom, pos, cn
    Optional: baf
    실제 파이프라인 출력 (긴 컬럼 포함) 자동 처리
    """
    df = _read_table(path)
    df.columns = [str(c).strip().lower() for c in df.columns]

    aliases = {
        'chromosome': 'chrom', 'chr': 'chrom',
        'position': 'pos', 'start': 'pos', 'chromstart': 'pos', 'chrom_start': 'pos',
        'copy_number': 'cn', 'copynumber': 'cn', 'copy-number': 'cn',
        'copy_number_signal': 'cn',                 # ← 파이프라인 실제 컬럼
        'log2_corrected': 'cn_log2',               # ← 참고용 보존
        'vaf': 'baf', 'b_allele_freq': 'baf',
        'baf_value': 'baf', 'allele_freq': 'baf',
        'bin_baf': 'baf',                           # ← 파이프라인 실제 컬럼
    }
    for src, dst in aliases.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})

    missing = [c for c in ['chrom', 'pos', 'cn'] if c not in df.columns]
    if missing:
        raise ValueError(f'signal file missing required columns: {missing}')

    df['chrom'] = _norm_chrom(df['chrom'])
    df['pos']   = pd.to_numeric(df['pos'],  errors='coerce')
    df['cn']    = pd.to_numeric(df['cn'],   errors='coerce')
    if 'baf' in df.columns:
        df['baf'] = pd.to_numeric(df['baf'], errors='coerce')

    keep = ['chrom', 'pos', 'cn'] + (['baf'] if 'baf' in df.columns else [])
    return df[keep].dropna(subset=['chrom', 'pos', 'cn']).sort_values(['chrom', 'pos'])


def load_cnv(path):
    df = _read_table(path)
    df.columns = [str(c).strip().lower() for c in df.columns]

    aliases = {
        'chromosome': 'chrom', 'chr': 'chrom',
        'stop': 'end', 'chromend': 'end', 'chrom_end': 'end',
        'chromstart': 'start', 'chrom_start': 'start',
        'copy_number': 'cn', 'copynumber': 'cn',
        'event': 'type', 'svtype': 'type',
        'label': 'name',
    }
    for src, dst in aliases.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})

    missing = [c for c in ['chrom', 'start', 'end'] if c not in df.columns]
    if missing:
        raise ValueError(f'cnv file missing required columns: {missing}')

    df['chrom'] = _norm_chrom(df['chrom'])
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end']   = pd.to_numeric(df['end'],   errors='coerce')
    if 'cn' in df.columns:
        df['cn'] = pd.to_numeric(df['cn'], errors='coerce')

    for col in ['type', 'name', 'iscn', 'color']:
        if col not in df.columns:
            df[col] = None

    return df.dropna(subset=['chrom', 'start', 'end']).sort_values(['chrom', 'start'])


def load_genes(path):
    if not path:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'color'])
    df = _read_table(path)
    df.columns = [str(c).strip().lower() for c in df.columns]
    aliases = {'chromosome': 'chrom', 'chr': 'chrom', 'stop': 'end', 'gene': 'name'}
    for src, dst in aliases.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})
    missing = [c for c in ['chrom', 'start', 'end', 'name'] if c not in df.columns]
    if missing:
        raise ValueError(f'gene file missing required columns: {missing}')
    df['chrom'] = _norm_chrom(df['chrom'])
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end']   = pd.to_numeric(df['end'],   errors='coerce')
    if 'color' not in df.columns:
        df['color'] = None
    return df[['chrom','start','end','name','color']].dropna(subset=['chrom','start','end','name'])
