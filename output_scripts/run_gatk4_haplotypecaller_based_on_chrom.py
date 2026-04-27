#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = HaplotypeCaller_based_on_chrom

"""
Tool: gatk4 (4.4.0.0)
Profile: HaplotypeCaller_based_on_chrom
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4 Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='Sequence identifier for sample tracking (Default: )')
    parser.add_argument('--Chromosome', required=True, default='', help='Target chromosome or interval (L option) (Default: )')
    parser.add_argument('--BamDir', required=True, default='', help='Directory containing the input BAM file (Default: )')
    parser.add_argument('--ResultDir', required=True, default='', help='Output directory for VCF and GVCF files (Default: )')
    parser.add_argument('--ReferenceFasta', required=False, default='/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa', help='Reference genome FASTA (Default: /storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa)')
    parser.add_argument('--DbsnpVcf', required=False, default='/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz', help='Known variation database (dbSNP) (Default: /storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz)')
    parser.add_argument('--InputSuffix', required=False, default='analysisReady', help='Suffix of input BAM (e.g., analysisReady, recal, dedup) (Default: analysisReady)')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='Path to singularity executable (Default: singularity)')
    parser.add_argument('--gatk_bin', required=False, default='gatk', help='GATK binary path inside container (Default: gatk)')
    parser.add_argument('--sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='GATK singularity image file (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='Mount paths for singularity (Default: /storage,/data)')
    parser.add_argument('--Threads', required=False, default='4', help='Threads for ParallelGC (Default: 4)')
    parser.add_argument('--Memory', required=False, default='32g', help='Java Xmx memory setting (e.g., 32g) (Default: 32g)')
    parser.add_argument('--Ploidy', required=False, default='2', help='Sample ploidy (Default: 2) (Default: 2)')
    parser.add_argument('--TmpDir', required=False, default='[ResultDir]/tmp/[SeqID]_[Chromosome]', help='Temporary directory for GATK (Default: [ResultDir]/tmp/[SeqID]_[Chromosome])')
    parser.add_argument('--IncludeNonVariant', required=False, default='false', help='Whether to include non-variant sites in GenotypeGVCFs (Default: false)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    Chromosome = args.Chromosome
    BamDir = args.BamDir
    ResultDir = args.ResultDir
    ReferenceFasta = args.ReferenceFasta
    DbsnpVcf = args.DbsnpVcf
    InputSuffix = args.InputSuffix
    singularity_bin = args.singularity_bin
    gatk_bin = args.gatk_bin
    sif = args.sif
    bind = args.bind
    Threads = args.Threads
    Memory = args.Memory
    Ploidy = args.Ploidy
    TmpDir = args.TmpDir
    IncludeNonVariant = args.IncludeNonVariant

    # --- [Output Paths] ---
    OutGvcf = f"{ResultDir}/{SeqID}.{Chromosome}.gvcf.gz"
    OutVcf = f"{ResultDir}/{SeqID}.{Chromosome}.vcf.gz"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {sif} {gatk_bin} --java-options '-XX:ParallelGCThreads={Threads} -Xmx{Memory} -Djava.io.tmpdir={TmpDir}'  HaplotypeCaller  -R {ReferenceFasta}  -I {BamDir}/{SeqID}.{InputSuffix}.bam  -L {Chromosome}  -ploidy {Ploidy}  -stand-call-conf 30  --dbsnp {DbsnpVcf}  -O {OutGvcf}  -ERC GVCF && 
{singularity_bin} exec -B {bind} {sif} {gatk_bin} --java-options '-XX:ParallelGCThreads={Threads} -Xmx{Memory} -Djava.io.tmpdir={TmpDir}'  GenotypeGVCFs  --include-non-variant-sites {IncludeNonVariant}  -R {ReferenceFasta}  -V {OutGvcf}  -O {OutVcf}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(BamDir) if '.' in os.path.basename(BamDir) else BamDir, exist_ok=True)
    os.makedirs(os.path.dirname(ResultDir) if '.' in os.path.basename(ResultDir) else ResultDir, exist_ok=True)
    os.makedirs(os.path.dirname(TmpDir) if '.' in os.path.basename(TmpDir) else TmpDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()