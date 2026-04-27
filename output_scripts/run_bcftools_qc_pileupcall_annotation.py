#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = bcftools
# VERSION = 1.23
# THREADS = 1
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
    InputSuffix = args.InputSuffix
    bcftools_bin = args.bcftools_bin
    Threads = args.Threads
    MinBQ = args.MinBQ
    MinMQ = args.MinMQ
    AnnotationQuery = args.AnnotationQuery

    # --- [Output Paths] ---
    raw_vcf = f"{ResultDir}/{SeqID}.{Chromosome}.{InputSuffix}.raw.vcf.gz"
    ann_vcf = f"{ResultDir}/{SeqID}.{Chromosome}.{InputSuffix}.ann.vcf.gz"

    # --- [Command Execution] ---
    cmd = f"{bcftools_bin} mpileup -f {ReferenceFasta} -T {SitesVcfGz} -r {Chromosome} -q {MinMQ} -Q {MinBQ} -a FORMAT/AD,FORMAT/DP -Ou --threads {Threads} {BamDir}/{SeqID}.{InputSuffix}.bam | {bcftools_bin} call -Am -Oz -o {raw_vcf} --threads {Threads} && {bcftools_bin} index -f {raw_vcf} --threads {Threads} && {bcftools_bin} annotate -a {PopAfAnnotVcf} -c {AnnotationQuery} -h {PopAfHeaderHdr} -Oz -o {ann_vcf} {raw_vcf} --threads {Threads} && {bcftools_bin} index -f {ann_vcf} --threads {Threads}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(ResultDir) if '.' in os.path.basename(ResultDir) else ResultDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()