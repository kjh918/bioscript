#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = Ginkgo_Final_Optimized
# VERSION = 1.0.0
# THREADS = 1
# PROFILE = absolute_path_runner

"""
Tool: Ginkgo_Final_Optimized (1.0.0)
Profile: absolute_path_runner
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Ginkgo_Final_Optimized Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--BinSize', required=True, default='', help='No description (Default: )')
    parser.add_argument('--ReadLength', required=True, default='', help='No description (Default: )')
    parser.add_argument('--Threads', required=False, default='1', help='No description (Default: 1)')
    parser.add_argument('--GinkgoHomeDir', required=False, default='/storage/apps/ginkgo', help='No description (Default: /storage/apps/ginkgo)')
    parser.add_argument('--BamToBedDir', required=False, default='/data/cbNIPT/bamToBeds', help='No description (Default: /data/cbNIPT/bamToBeds)')
    parser.add_argument('--NormalDir', required=False, default='/data/cbNIPT/bamToBeds', help='No description (Default: /data/cbNIPT/bamToBeds)')
    parser.add_argument('--ResultBaseDir', required=False, default='/data/cbNIPT/ginkgo_analysis', help='No description (Default: /data/cbNIPT/ginkgo_analysis)')
    parser.add_argument('--Genome', required=False, default='hg38', help='No description (Default: hg38)')
    parser.add_argument('--NormalSamples', required=False, default='Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz', help='No description (Default: Normal_01.bed.gz Normal_02.bed.gz Normal_03.bed.gz)')
    parser.add_argument('--R_Args_Base', required=False, default='status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1', help='No description (Default: status.xml data 2 [BinMeth] ward euclidean 1 refDummy.bed_mapped 0 ploidyDummy.txt 0 1)')
    parser.add_argument('--R_Args_Reclust', required=False, default='status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0', help='No description (Default: status.xml [BinMeth] ward euclidean 0 ploidyDummy.txt 0)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    BinSize = args.BinSize
    ReadLength = args.ReadLength
    Threads = args.Threads
    GinkgoHomeDir = args.GinkgoHomeDir
    BamToBedDir = args.BamToBedDir
    NormalDir = args.NormalDir
    ResultBaseDir = args.ResultBaseDir
    Genome = args.Genome
    NormalSamples = args.NormalSamples
    R_Args_Base = args.R_Args_Base
    R_Args_Reclust = args.R_Args_Reclust

    # --- [Output Paths] ---
    BinMeth = f"variable_{BinSize}000_{ReadLength}_bwa"
    WorkDir = f"{ResultBaseDir}/{SeqID}/{SeqID}_{BinSize}kb"
    BinUnsorted = f"{GinkgoHomeDir}/scripts/binUnsorted"
    BinFile = f"{GinkgoHomeDir}/genomes/{Genome}/original/{BinMeth}"
    BinCount = f"$(wc -l < {BinFile})"
    ProcessR = f"{GinkgoHomeDir}/scripts/process.R"
    ReclustR = f"{GinkgoHomeDir}/scripts/reclust.R"
    CNVCaller = f"{GinkgoHomeDir}/scripts/CNVcaller"

    # --- [Command Execution] ---
    cmd = f"mkdir -p {WorkDir} && ln -sf {NormalDir}/*.bed.gz {WorkDir}/ && ln -sf {BamToBedDir}/{SeqID}.bed.gz {WorkDir}/ &&
find {WorkDir} -maxdepth 1 -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {BinUnsorted} {BinFile} {BinCount} <(zcat '{}' | awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}') $N '{}'_mapped' &&
find {WorkDir} -maxdepth 1 -name '*.bed' ! -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {BinUnsorted} {BinFile} {BinCount} <(awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}' '{}') $N '{}'_mapped' &&
paste {WorkDir}/*_mapped > {WorkDir}/data &&
{ProcessR} {GinkgoHomeDir}/genomes/{Genome}/original {WorkDir} {R_Args_Base} && {ReclustR} {GinkgoHomeDir}/genomes/{Genome}/original {WorkDir} {R_Args_Reclust} && {CNVCaller} {WorkDir}/SegCopy {WorkDir}/CNV1 {WorkDir}/CNV2"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(GinkgoHomeDir) if '.' in os.path.basename(GinkgoHomeDir) else GinkgoHomeDir, exist_ok=True)
    os.makedirs(os.path.dirname(BamToBedDir) if '.' in os.path.basename(BamToBedDir) else BamToBedDir, exist_ok=True)
    os.makedirs(os.path.dirname(NormalDir) if '.' in os.path.basename(NormalDir) else NormalDir, exist_ok=True)
    os.makedirs(os.path.dirname(ResultBaseDir) if '.' in os.path.basename(ResultBaseDir) else ResultBaseDir, exist_ok=True)
    os.makedirs(os.path.dirname(WorkDir) if '.' in os.path.basename(WorkDir) else WorkDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()