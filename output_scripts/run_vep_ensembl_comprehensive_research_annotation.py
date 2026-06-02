#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = vep_ensembl
# VERSION = 110
# THREADS = 8
# PROFILE = comprehensive_research_annotation

"""
Tool: vep_ensembl (110)
Profile: comprehensive_research_annotation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="vep_ensembl Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--VcfTag', required=True, default='', help='Input VCF tag (e.g. mutect2.filtered) (Default: )')
    parser.add_argument('--vepCacheDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--genomeFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vepBam', required=True, default='', help='BAM file for transcript-matching (Default: )')
    parser.add_argument('--clinvarData', required=True, default='', help='No description (Default: )')
    parser.add_argument('--cosmicData', required=False, default='', help='No description (Default: )')
    parser.add_argument('--alphaMissense', required=True, default='', help='No description (Default: )')
    parser.add_argument('--caddSnp', required=True, default='', help='No description (Default: )')
    parser.add_argument('--caddIndel', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vep_vcf', required=False, default='[vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf', help='No description (Default: [vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf)')
    parser.add_argument('--vep_stats', required=False, default='[vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf_summary.txt', help='No description (Default: [vcfDir]/[SeqID].[VcfTag].vep.ensembl.vcf_summary.txt)')
    parser.add_argument('--assembly', required=False, default='GRCh38', help='GRCh38 or GRCh37 (Default: GRCh38)')
    parser.add_argument('--species', required=False, default='homo_sapiens', help='No description (Default: homo_sapiens)')
    parser.add_argument('--Threads', required=False, default='8', help='No description (Default: 8)')
    parser.add_argument('--buffer_size', required=False, default='50000', help='No description (Default: 50000)')
    parser.add_argument('--cosmic_args', required=False, default='', help='No description (Default: )')
    parser.add_argument('--clinvar_args', required=False, default='--custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN', help='No description (Default: --custom file=[clinvarData],short_name=ClinVar,format=vcf,fields=CLNSIG%CLNDN%ORIGIN)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--vep_sif', required=False, default='/storage/images/vep-110.sif', help='No description (Default: /storage/images/vep-110.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    vcfDir = args.vcfDir
    VcfTag = args.VcfTag
    vepCacheDir = args.vepCacheDir
    genomeFasta = args.genomeFasta
    vepBam = args.vepBam
    clinvarData = args.clinvarData
    cosmicData = args.cosmicData
    alphaMissense = args.alphaMissense
    caddSnp = args.caddSnp
    caddIndel = args.caddIndel
    vep_vcf = args.vep_vcf
    vep_stats = args.vep_stats
    assembly = args.assembly
    species = args.species
    Threads = args.Threads
    buffer_size = args.buffer_size
    cosmic_args = args.cosmic_args
    clinvar_args = args.clinvar_args
    singularity_bin = args.singularity_bin
    vep_sif = args.vep_sif
    bind = args.bind

    # --- [Output Paths] ---
    if not vep_vcf:
    vep_vcf = f"{vcfDir}/{SeqID}.{VcfTag}.vep.ensembl.vcf"
    if not vep_stats:
    vep_stats = f"{vcfDir}/{SeqID}.{VcfTag}.vep.ensembl.vcf_summary.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {vep_sif} vep --force_overwrite --offline --cache --dir_cache {vepCacheDir} --fasta {genomeFasta} --bam {vepBam} --species {species} --assembly {assembly} --input_file {vcfDir}/{SeqID}.{VcfTag}.vcf --output_file {vep_vcf} --vcf --stats_text --hgvs --hgvsg --canonical --exclude_predicted --ccds --uniprot --domains --fork {Threads} --buffer_size {buffer_size} {clinvar_args} {cosmic_args} --plugin AlphaMissense,file={alphaMissense} --plugin CADD,snv={caddSnp},indels={caddIndel}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if vcfDir:
        _tgt = os.path.dirname(vcfDir) if os.path.splitext(vcfDir)[1] else vcfDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if vep_vcf:
        _tgt = os.path.dirname(vep_vcf) if os.path.splitext(vep_vcf)[1] else vep_vcf
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if vep_stats:
        _tgt = os.path.dirname(vep_stats) if os.path.splitext(vep_stats)[1] else vep_stats
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if vepCacheDir:
        _tgt = os.path.dirname(vepCacheDir) if os.path.splitext(vepCacheDir)[1] else vepCacheDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()