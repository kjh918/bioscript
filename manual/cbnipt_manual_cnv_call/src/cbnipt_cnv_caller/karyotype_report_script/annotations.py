import math
import pandas as pd
from config import DEFAULT_COLORS


def infer_event_type(row: pd.Series) -> str:
    raw = str(row.get('type') or '').lower()
    name = str(row.get('name') or row.get('iscn') or '').lower()
    cn = row.get('cn')

    if any(k in raw + ' ' + name for k in ['loss', 'del', 'mono']):
        return 'loss'
    if any(k in raw + ' ' + name for k in ['gain', 'dup', 'tri', 'amp']):
        return 'gain'
    try:
        cn = float(cn)
        if cn < 2:
            return 'loss'
        if cn > 2:
            return 'gain'
    except Exception:
        pass
    return 'neutral'


def cnv_label(row: pd.Series) -> str:
    for key in ['name', 'iscn']:
        v = row.get(key)
        if v is not None and str(v).strip() and str(v).lower() != 'nan':
            return str(v)
    kind = infer_event_type(row)
    cn = row.get('cn')
    cn_txt = ''
    try:
        if not math.isnan(float(cn)):
            cn_txt = f' CN={float(cn):g}'
    except Exception:
        pass
    return f'{kind.upper()}{cn_txt}'


def build_ideogram_annotations(cnv: pd.DataFrame, genes: pd.DataFrame) -> list[dict]:
    """Build annotation objects accepted directly by Ideogram.js.

    CNV = rectangle, Gene = circle.
    Coordinates are converted to 0-based starts for Ideogram.js display.
    """
    annots: list[dict] = []
    for _, row in cnv.iterrows():
        kind = infer_event_type(row)
        color = row.get('color')
        if not color or str(color).lower() == 'nan':
            color = DEFAULT_COLORS[kind]
        annots.append({
            'name': cnv_label(row),
            'chr': str(row['chrom']),
            'start': max(0, int(row['start']) - 1),
            'stop': max(0, int(row['end']) - 1),
            'color': color,
            'shape': 'rectangle',
            'trackIndex': 0,
        })

    for _, row in genes.iterrows():
        color = row.get('color')
        if not color or str(color).lower() == 'nan':
            color = DEFAULT_COLORS['gene']
        annots.append({
            'name': str(row['name']),
            'chr': str(row['chrom']),
            'start': max(0, int(row['start']) - 1),
            'stop': max(0, int(row['end']) - 1),
            'color': color,
            'shape': 'circle',
            'trackIndex': 1,
        })
    return annots


def summarize_cnv(cnv: pd.DataFrame) -> list[dict]:
    rows = []
    for _, row in cnv.iterrows():
        rows.append({
            'chrom': str(row['chrom']),
            'start': int(row['start']),
            'end': int(row['end']),
            'length_mb': max(0, int(row['end']) - int(row['start'])) / 1e6,
            'type': infer_event_type(row),
            'cn': None if pd.isna(row.get('cn')) else row.get('cn'),
            'name': cnv_label(row),
        })
    return rows
