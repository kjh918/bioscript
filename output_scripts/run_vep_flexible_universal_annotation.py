#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = vep_flexible
# VERSION = 110
# THREADS = 1
# PROFILE = universal_annotation

"""
Tool: vep_flexible (110)
Profile: universal_annotation
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="vep_flexible Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vcfDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--VcfTag', required=True, default='', help='Input VCF tag (e.g. mutect2.bias.filtered) (Default: )')
    parser.add_argument('--vepCacheDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--vepPluginDir', required=True, default='', help='No description (Default: )')
    parser.add_argument('--genomeFasta', required=True, default='', help='No description (Default: )')
    parser.add_argument('--clinvarData', required=True, default='', help='ClinVar VCF path (Assembly-specific) (Default: )')
    parser.add_argument('--cosmicData', required=False, default='', help='COSMIC VCF path (Assembly-specific, optional) (Default: )')
    parser.add_argument('--alphaMissense', required=True, default='', help='No description (Default: )')
    parser.add_argument('--caddSnp', required=True, default='', help='No description (Default: )')
    parser.add_argument('--caddIndel', required=True, default='', help='No description (Default: )')
    parser.add_argument('--assembly', required=False, default='GRCh38', help='GRCh38 or GRCh37 (Default: GRCh38)')
    parser.add_argument('--species', required=False, default='homo_sapiens', help='No description (Default: homo_sapiens)')
    parser.add_argument('--Threads', required=False, default='8', help='No description (Default: 8)')
    parser.add_argument('--buffer_size', required=False, default='50000', help='No description (Default: 50000)')
    parser.add_argument('--cosmic_args', required=False, default='', help='Custom flag for COSMIC (e.g. --custom file=...) (Default: )')
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
    vepPluginDir = args.vepPluginDir
    genomeFasta = args.genomeFasta
    clinvarData = args.clinvarData
    cosmicData = args.cosmicData
    alphaMissense = args.alphaMissense
    caddSnp = args.caddSnp
    caddIndel = args.caddIndel
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
    vep_vcf = f"{vcfDir}/{SeqID}.{VcfTag}.vep.refseq.vcf"
    vep_stats = f"{vcfDir}/{SeqID}.{VcfTag}.vep.refseq.vcf_summary.txt"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {vep_sif} vep --force_overwrite --offline --cache --dir_cache {vepCacheDir} --dir_plugins {vepPluginDir} --fasta {genomeFasta} --species {species} --assembly {assembly} --input_file {vcfDir}/{SeqID}.{VcfTag}.vcf --output_file {vep_vcf} --vcf --stats_text --refseq --show_ref_allele --uploaded_allele --use_transcript_ref --variant_class --sift b --polyphen b --gene_phenotype --numbers --hgvs --hgvsg --symbol --canonical --biotype --regulatory --mirna --check_existing --max_af --af_1kg --af_gnomade --exclude_predicted --pick --flag_pick --fork {Threads} --buffer_size {buffer_size} {clinvar_args} {cosmic_args} --plugin AlphaMissense,file={alphaMissense} --plugin CADD,snv={caddSnp},indels={caddIndel}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(vcfDir) if '.' in os.path.basename(vcfDir) else vcfDir, exist_ok=True)
    os.makedirs(os.path.dirname(vepCacheDir) if '.' in os.path.basename(vepCacheDir) else vepCacheDir, exist_ok=True)
    os.makedirs(os.path.dirname(vepPluginDir) if '.' in os.path.basename(vepPluginDir) else vepPluginDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()