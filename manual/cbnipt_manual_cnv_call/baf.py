import pandas as pd
import numpy as np
import pysam

def count_baf_for_sites(bam_path, sites_df, min_mapq, min_baseq):
    bam = pysam.AlignmentFile(bam_path, "rb")
    res = []
    for row in sites_df.itertuples():
        rc, ac = 0, 0
        for col in bam.pileup(row.chrom, row.pos-1, row.pos, truncate=True, min_base_quality=min_baseq):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip: continue
                rd = pr.alignment
                if rd.mapping_quality < min_mapq or rd.is_duplicate: continue
                base = rd.query_sequence[pr.query_position].upper()
                if base == str(row.ref).upper(): rc += 1
                elif base == str(row.alt).upper(): ac += 1
        d = rc + ac
        res.append({**row._asdict(), "ref_count": rc, "alt_count": ac, "depth": d, "baf": ac/d if d>0 else np.nan})
    bam.close()
    return pd.DataFrame(res)