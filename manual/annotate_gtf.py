#!/usr/bin/env python3

import os
import re
import argparse
import pysam
import pandas as pd


def strip_ensembl_version(ensembl_id: str) -> str:
    return re.sub(r"\.\d+$", "", str(ensembl_id))


def parse_gtf_attributes(attr_string: str) -> dict:
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item:
            continue

        m = re.match(r'(\S+)\s+"([^"]+)"', item)
        if m:
            attrs[m.group(1)] = m.group(2)

    return attrs


def detect_ensembl_type(clean_id: str) -> str:
    if clean_id.startswith("ENSG"):
        return "gene"
    if clean_id.startswith("ENST"):
        return "transcript"
    return "unknown"


def build_gtf_index_from_tabix(gtf_gz: str) -> dict:
    if not os.path.exists(gtf_gz):
        raise FileNotFoundError(gtf_gz)

    if not os.path.exists(gtf_gz + ".tbi"):
        raise FileNotFoundError(
            f"Tabix index not found: {gtf_gz}.tbi\n"
            f"Create index first: tabix -p gff {gtf_gz}"
        )

    index = {
        "gene": {},
        "transcript": {}
    }

    with pysam.TabixFile(gtf_gz) as tbx:
        for contig in tbx.contigs:
            for line in tbx.fetch(contig):
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attrs_raw = fields

                if feature not in ["gene", "transcript"]:
                    continue

                attrs = parse_gtf_attributes(attrs_raw)

                gene_id = attrs.get("gene_id", "")
                transcript_id = attrs.get("transcript_id", "")

                gene_id_clean = strip_ensembl_version(gene_id)
                transcript_id_clean = strip_ensembl_version(transcript_id)

                record = {
                    "gtf_feature": feature,
                    "ensembl_gene_id": gene_id_clean,
                    "ensembl_gene_id_version": gene_id,
                    "ensembl_transcript_id": transcript_id_clean,
                    "ensembl_transcript_id_version": transcript_id,
                    "gene_symbol": attrs.get("gene_name", ""),
                    "gene_biotype": attrs.get("gene_biotype", attrs.get("gene_type", "")),
                    "transcript_biotype": attrs.get("transcript_biotype", attrs.get("transcript_type", "")),
                    "chromosome": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "source": source
                }

                if feature == "gene" and gene_id_clean:
                    index["gene"][gene_id_clean] = record

                if feature == "transcript" and transcript_id_clean:
                    index["transcript"][transcript_id_clean] = record

    print(f"[INFO] Loaded genes: {len(index['gene'])}")
    print(f"[INFO] Loaded transcripts: {len(index['transcript'])}")

    return index


def get_annotation(index: dict, ensembl_id: str) -> dict:
    clean_id = strip_ensembl_version(ensembl_id)
    id_type = detect_ensembl_type(clean_id)

    if id_type == "gene":
        info = index["gene"].get(clean_id)
    elif id_type == "transcript":
        info = index["transcript"].get(clean_id)
    else:
        info = index["gene"].get(clean_id) or index["transcript"].get(clean_id)

    if info is None:
        return {
            "clean_ensembl_id": clean_id,
            "ensembl_id_type": id_type,
            "matched": False,
            "gtf_feature": "",
            "ensembl_gene_id": "",
            "ensembl_gene_id_version": "",
            "ensembl_transcript_id": "",
            "ensembl_transcript_id_version": "",
            "gene_symbol": "",
            "gene_biotype": "",
            "transcript_biotype": "",
            "chromosome": "",
            "start": "",
            "end": "",
            "strand": "",
            "source": ""
        }

    out = {
        "clean_ensembl_id": clean_id,
        "ensembl_id_type": id_type,
        "matched": True
    }
    out.update(info)
    return out


def read_table_auto(input_file: str) -> pd.DataFrame:
    if input_file.endswith(".csv"):
        return pd.read_csv(input_file, dtype=str)
    return pd.read_csv(input_file, sep="\t", dtype=str)


def annotate_input_table(
    input_file: str,
    id_col: str,
    gtf_gz: str,
    output_tsv: str
) -> pd.DataFrame:

    df = read_table_auto(input_file)

    if id_col not in df.columns:
        raise ValueError(
            f"ID column not found: {id_col}\n"
            f"Available columns: {', '.join(df.columns)}"
        )

    index = build_gtf_index_from_tabix(gtf_gz)

    annotation_records = []

    for x in df[id_col].fillna("").astype(str):
        annotation_records.append(get_annotation(index, x))

    anno_df = pd.DataFrame(annotation_records)

    merged = pd.concat(
        [
            df.reset_index(drop=True),
            anno_df.reset_index(drop=True)
        ],
        axis=1
    )

    front_cols = [
        id_col,
        "clean_ensembl_id",
        "ensembl_id_type",
        "matched",
        "gene_symbol",
        "ensembl_gene_id",
        "ensembl_gene_id_version",
        "ensembl_transcript_id",
        "ensembl_transcript_id_version",
        "gene_biotype",
        "transcript_biotype",
        "gtf_feature",
        "chromosome",
        "start",
        "end",
        "strand",
        "source"
    ]

    existing_front_cols = [c for c in front_cols if c in merged.columns]
    other_cols = [c for c in merged.columns if c not in existing_front_cols]

    merged = merged[existing_front_cols + other_cols]

    merged.to_csv(output_tsv, sep="\t", index=False)

    print(f"[INFO] Input rows: {len(df)}")
    print(f"[INFO] Matched rows: {merged['matched'].sum()}")
    print(f"[INFO] Output saved: {output_tsv}")

    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Merge ENSEMBL gene/transcript annotation from bgzip+tabix indexed GTF into input TSV/CSV"
    )

    parser.add_argument(
        "-g", "--gtf",
        required=True,
        help="bgzip-compressed and tabix-indexed GTF file"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV/CSV file"
    )

    parser.add_argument(
        "-c", "--id-col",
        required=True,
        help="Column name containing ENSEMBL gene/transcript IDs"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output merged TSV file"
    )

    args = parser.parse_args()

    annotate_input_table(
        input_file=args.input,
        id_col=args.id_col,
        gtf_gz=args.gtf,
        output_tsv=args.output
    )


if __name__ == "__main__":
    main()