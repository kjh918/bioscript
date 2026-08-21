"""
BAM fragment 길이 분포 확인 스크립트
python check_bam_fragments.py --bam sample.bam --chrom chr1 --n 100000
"""
import argparse
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam",   required=True)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--start", type=int, default=1000000)
    parser.add_argument("--end",   type=int, default=2000000)
    parser.add_argument("--n",     type=int, default=100000)
    parser.add_argument("--out",   default="bam_check.png")
    args = parser.parse_args()

    lengths = []
    seen = set()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(args.chrom, args.start, args.end):
            if len(lengths) >= args.n:
                break
            if (read.is_unmapped or read.is_duplicate or
                    read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < 20:
                continue
            if read.is_paired and not read.is_read1:
                continue
            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)

            if read.is_paired and read.template_length != 0:
                flen = abs(read.template_length)
            else:
                flen = read.query_length or 0
            if flen > 0:
                lengths.append(flen)

    lengths = np.array(lengths)
    print(f"=== Fragment 길이 분포 (n={len(lengths):,}) ===")
    print(f"  median : {np.median(lengths):.0f} bp")
    print(f"  mean   : {np.mean(lengths):.0f} bp")
    print(f"  range  : [{lengths.min()}, {lengths.max()}]")
    for thr in [35, 80, 120, 150, 180, 250, 500]:
        pct = (lengths <= thr).mean() * 100
        print(f"  <= {thr:3d} bp : {pct:.1f}%  ({(lengths <= thr).sum():,})")

    print(f"\nWPS 범위별:")
    print(f"  S-WPS (35-80bp)  : {((lengths>=35)&(lengths<=80)).sum():,}  ({((lengths>=35)&(lengths<=80)).mean()*100:.1f}%)")
    print(f"  L-WPS (120-180bp): {((lengths>=120)&(lengths<=180)).sum():,}  ({((lengths>=120)&(lengths<=180)).mean()*100:.1f}%)")
    print(f"  SHORT (<=150bp)  : {(lengths<=150).sum():,}  ({(lengths<=150).mean()*100:.1f}%)")
    print(f"  LONG  (>=151bp)  : {(lengths>=151).sum():,}  ({(lengths>=151).mean()*100:.1f}%)")

    fig, ax = plt.subplots(figsize=(12, 5), facecolor="white")
    ax.hist(lengths, bins=200, range=(0, 600), color="#3278d6", alpha=0.7)
    for x, c in [(35,"g"),(80,"g"),(120,"r"),(150,"k"),(180,"r")]:
        ax.axvline(x, color=c, lw=1.2, ls="--", label=f"{x}bp")
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Count")
    ax.set_title(f"Fragment length distribution  {args.chrom}:{args.start}-{args.end}")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"\n그래프 저장: {args.out}")

if __name__ == "__main__":
    main()