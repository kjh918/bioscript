#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = vcf2maf
# VERSION = 1.6.21
# THREADS = 1
# PROFILE = maf_conversion

"""
Tool: vcf2maf (1.6.21)
Profile: maf_conversion
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="vcf2maf Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--VcfTag', required=True, default='', help='Input VCF tag (e.g. mutect2.filtered.vep.refseq) (Default: )')
    parser.add_argument('--NormalID', required=True, default='NORMAL', help='Normal sample ID (Default: NORMAL)')
    parser.add_argument('--genomeFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--retain_info', required=False, default='DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus', help='MAF에 컬럼으로 남길 VCF INFO 필드 리스트 (Default: DP,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,PON,ROQ,TLOD,FILTER,pArtifact,artiStatus)')
    parser.add_argument('--retain_ann', required=False, default='HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER', help='MAF에 컬럼으로 남길 VEP ANN(CSQ) 필드 리스트 (Default: HGVSg,am_class,am_pathogenicity,CADD_PHRED,CADD_RAW,ClinVar_CLNSIG,ClinVar_CLNDN,ClinVar_ORIGIN,gnomADe_AF,gnomADe_AFR_AF,gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,MAX_AF,MAX_AF_POPS,TRANSCRIPTION_FACTORS,COSMIC,COSMIC_HGVSC,COSMIC_LEGACY_ID,COSMIC_TIER)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--vcf2maf_sif', required=False, default='/storage/images/vcf2maf.sif', help='No description (Default: /storage/images/vcf2maf.sif)')
    parser.add_argument('--vcf2maf_bin', required=False, default='vcf2maf.pl', help='No description (Default: vcf2maf.pl)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    vcfDir = args.vcfDir
    VcfTag = args.VcfTag
    NormalID = args.NormalID
    genomeFasta = args.genomeFasta
    retain_info = args.retain_info
    retain_ann = args.retain_ann
    singularity_bin = args.singularity_bin
    vcf2maf_sif = args.vcf2maf_sif
    vcf2maf_bin = args.vcf2maf_bin
    bind = args.bind

    # --- [Output Paths] ---
    output_maf = f"{vcfDir}/{SeqID}.{VcfTag}.maf"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {vcf2maf_sif} {vcf2maf_bin} --inhibit-vep --input-vcf {vcfDir}/{SeqID}.{VcfTag}.vcf --output-maf {output_maf} --tumor-id {SeqID} --normal-id {NormalID} --ref-fasta {genomeFasta} --retain-info {retain_info} --retain-ann {retain_ann}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(vcfDir) if '.' in os.path.basename(vcfDir) else vcfDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()