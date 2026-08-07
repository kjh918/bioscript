#!/usr/bin/env python3
import argparse, sys
from pathlib import Path

from io_utils import load_signal, load_cnv, load_genes
from annotations import build_ideogram_annotations, summarize_cnv
from plots import build_plot_html
from html_report import render_report


def parse_args():
    p = argparse.ArgumentParser(description='Generate static Karyotype/CNV HTML report')
    p.add_argument('--signal',    required=True)
    p.add_argument('--cnv',       required=True)
    p.add_argument('--genes',     default=None)
    p.add_argument('--sample-id', default='SAMPLE')
    p.add_argument('--sex',       choices=['female','male'], default='female')
    p.add_argument('--output',    required=True)
    return p.parse_args()


def main():
    args   = parse_args()
    signal = load_signal(args.signal)
    cnv    = load_cnv(args.cnv)
    genes  = load_genes(args.genes)

    print(f'  signal : {len(signal):,} rows  chroms={sorted(signal["chrom"].unique())}')
    print(f'  cnv    : {len(cnv)} segments')
    print(f'  genes  : {len(genes)} entries')

    anns   = build_ideogram_annotations(cnv, genes)
    rows   = summarize_cnv(cnv)
    plots  = build_plot_html(signal, cnv, include_plotlyjs=True)

    out = render_report(
        output=args.output, sample_id=args.sample_id, sex=args.sex,
        annotations=anns, cnv_rows=rows, plot_html=plots,
    )
    print(f'✅  Report: {Path(out).resolve()}')


if __name__ == '__main__':
    main()
