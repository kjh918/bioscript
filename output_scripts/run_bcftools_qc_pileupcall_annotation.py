#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bcftools
# VERSION = 1.23
# THREADS = 4
# PROFILE = QC_PileupCall_Annotation

"""
Tool: bcftools (1.23)
Profile: QC_PileupCall_Annotation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="bcftools Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier for the sample (Default: )')
    parser.add_argument('--Chromosome', required=True, default='', help='Target chromosome or region (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM file (Default: )')
    parser.add_argument('--ResultDir', required=True, default='', help='Output directory for VCF results (Default: )')
    parser.add_argument('--ReferenceFasta', required=True, default='', help='Path to reference genome FASTA (Default: )')
    parser.add_argument('--SitesVcfGz', required=True, default='', help='VCF file containing specific sites to genotype (Default: )')
    parser.add_argument('--PopAfAnnotVcf', required=True, default='', help='VCF/Tabix file for population AF annotation (Default: )')
    parser.add_argument('--PopAfHeaderHdr', required=True, default='', help='Header file describing the annotation columns (Default: )')
    parser.add_argument('--raw_vcf', required=False, default='[ResultDir]/[SeqID].[Chromosome].[InputSuffix].raw.vcf.gz', help='Raw VCF file after pileup and calling (Default: [ResultDir]/[SeqID].[Chromosome].[InputSuffix].raw.vcf.gz)')
    parser.add_argument('--ann_vcf', required=False, default='[ResultDir]/[SeqID].[Chromosome].[InputSuffix].ann.vcf.gz', help='Annotated VCF file with population AF (Default: [ResultDir]/[SeqID].[Chromosome].[InputSuffix].ann.vcf.gz)')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of input BAM (e.g., analysisReady, recal) (Default: analysisReady)')
    parser.add_argument('--bcftools_bin', required=False, default='/storage/home/jhkim/Apps/bcftools/bcftools', help='Path to bcftools binary (Default: /storage/home/jhkim/Apps/bcftools/bcftools)')
    parser.add_argument('--Threads', required=False, default='4', help='Number of threads for compression/processing (Default: 4)')
    parser.add_argument('--MinBQ', required=False, default='20', help='Minimum base quality for mpileup (Default: 20)')
    parser.add_argument('--MinMQ', required=False, default='30', help='Minimum mapping quality for mpileup (Default: 30)')
    parser.add_argument('--AnnotationQuery', required=False, default='CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN', help='Columns to extract from annotation VCF (Default: CHROM,POS,REF,ALT,INFO/KOVA_AF,INFO/KOVA_AN,INFO/GNOMAD_AF,INFO/GNOMAD_AN)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    Chromosome = args.Chromosome
    BamDir = args.BamDir
    ResultDir = args.ResultDir
    ReferenceFasta = args.ReferenceFasta
    SitesVcfGz = args.SitesVcfGz
    PopAfAnnotVcf = args.PopAfAnnotVcf
    PopAfHeaderHdr = args.PopAfHeaderHdr
    raw_vcf = args.raw_vcf
    ann_vcf = args.ann_vcf
    InputSuffix = args.InputSuffix
    bcftools_bin = args.bcftools_bin
    Threads = args.Threads
    MinBQ = args.MinBQ
    MinMQ = args.MinMQ
    AnnotationQuery = args.AnnotationQuery

    # --- [Output Paths] ---
    if not raw_vcf:
    raw_vcf = f"{ResultDir}/{SeqID}.{Chromosome}.{InputSuffix}.raw.vcf.gz"
    if not ann_vcf:
    ann_vcf = f"{ResultDir}/{SeqID}.{Chromosome}.{InputSuffix}.ann.vcf.gz"

    # --- [Command Execution] ---
    cmd = f"{bcftools_bin} mpileup -f {ReferenceFasta} -T {SitesVcfGz} -r {Chromosome} -q {MinMQ} -Q {MinBQ} -a FORMAT/AD,FORMAT/DP -Ou --threads {Threads} {BamDir}/{SeqID}.{InputSuffix}.bam | {bcftools_bin} call -Am -Oz -o {raw_vcf} --threads {Threads} && {bcftools_bin} index -f {raw_vcf} --threads {Threads} && {bcftools_bin} annotate -a {PopAfAnnotVcf} -c {AnnotationQuery} -h {PopAfHeaderHdr} -Oz -o {ann_vcf} {raw_vcf} --threads {Threads} && {bcftools_bin} index -f {ann_vcf} --threads {Threads}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if raw_vcf:
        _tgt = os.path.dirname(raw_vcf) if os.path.splitext(raw_vcf)[1] else raw_vcf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamDir:
        _tgt = os.path.dirname(BamDir) if os.path.splitext(BamDir)[1] else BamDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if ResultDir:
        _tgt = os.path.dirname(ResultDir) if os.path.splitext(ResultDir)[1] else ResultDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if ann_vcf:
        _tgt = os.path.dirname(ann_vcf) if os.path.splitext(ann_vcf)[1] else ann_vcf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()