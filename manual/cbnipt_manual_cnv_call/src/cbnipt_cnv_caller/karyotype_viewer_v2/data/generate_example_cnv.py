"""
Generate example CN/BAF TSV files for each chromosome.
Run once to create data/example_cnv_chr21.tsv, etc.

Usage: python data/generate_example_cnv.py
"""
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

import pandas as pd
from karyotype_viewer.demo_sample import make_demo_sample
from karyotype_viewer.core.data import demo_dataframe

sample = make_demo_sample()
out = pathlib.Path(__file__).parent

for chrom in ["5", "17", "21", "X"]:
    df = demo_dataframe(sample, chrom)
    path = out / f"example_cnv_chr{chrom}.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"Written: {path}  ({len(df)} rows)")
