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
    parser.add_argument('--BinMeth', required=False, default='variable_[BinSize]000_[ReadLength]_bwa', help='No description (Default: variable_[BinSize]000_[ReadLength]_bwa)')
    parser.add_argument('--WorkDir', required=False, default='[ResultBaseDir]/[SeqID]/[SeqID]_[BinSize]kb', help='No description (Default: [ResultBaseDir]/[SeqID]/[SeqID]_[BinSize]kb)')
    parser.add_argument('--BinUnsorted', required=False, default='[GinkgoHomeDir]/scripts/binUnsorted', help='No description (Default: [GinkgoHomeDir]/scripts/binUnsorted)')
    parser.add_argument('--BinFile', required=False, default='[GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth]', help='No description (Default: [GinkgoHomeDir]/genomes/[Genome]/original/[BinMeth])')
    parser.add_argument('--BinCount', required=False, default='$(wc -l < [BinFile])', help='No description (Default: $(wc -l < [BinFile]))')
    parser.add_argument('--ProcessR', required=False, default='[GinkgoHomeDir]/scripts/process.R', help='No description (Default: [GinkgoHomeDir]/scripts/process.R)')
    parser.add_argument('--ReclustR', required=False, default='[GinkgoHomeDir]/scripts/reclust.R', help='No description (Default: [GinkgoHomeDir]/scripts/reclust.R)')
    parser.add_argument('--CNVCaller', required=False, default='[GinkgoHomeDir]/scripts/CNVcaller', help='No description (Default: [GinkgoHomeDir]/scripts/CNVcaller)')
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
    BinMeth = args.BinMeth
    WorkDir = args.WorkDir
    BinUnsorted = args.BinUnsorted
    BinFile = args.BinFile
    BinCount = args.BinCount
    ProcessR = args.ProcessR
    ReclustR = args.ReclustR
    CNVCaller = args.CNVCaller
    NormalSamples = args.NormalSamples
    R_Args_Base = args.R_Args_Base
    R_Args_Reclust = args.R_Args_Reclust

    # --- [Output Paths] ---
    if not BinMeth:
    BinMeth = f"variable_{BinSize}000_{ReadLength}_bwa"
    if not WorkDir:
    WorkDir = f"{ResultBaseDir}/{SeqID}/{SeqID}_{BinSize}kb"
    if not BinUnsorted:
    BinUnsorted = f"{GinkgoHomeDir}/scripts/binUnsorted"
    if not BinFile:
    BinFile = f"{GinkgoHomeDir}/genomes/{Genome}/original/{BinMeth}"
    if not BinCount:
    BinCount = f"$(wc -l < {BinFile})"
    if not ProcessR:
    ProcessR = f"{GinkgoHomeDir}/scripts/process.R"
    if not ReclustR:
    ReclustR = f"{GinkgoHomeDir}/scripts/reclust.R"
    if not CNVCaller:
    CNVCaller = f"{GinkgoHomeDir}/scripts/CNVcaller"

    # --- [Command Execution] ---
    cmd = f"mkdir -p {WorkDir} && ln -sf {NormalDir}/*.bed.gz {WorkDir}/ && ln -sf {BamToBedDir}/{SeqID}.bed.gz {WorkDir}/ &&
find {WorkDir} -maxdepth 1 -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {BinUnsorted} {BinFile} {BinCount} <(zcat '{}' | awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}') $N '{}'_mapped' &&
find {WorkDir} -maxdepth 1 -name '*.bed' ! -name '*.gz' | xargs -I {} bash -c 'N=$(basename '{}' | sed 's/\.bed.*//'); {BinUnsorted} {BinFile} {BinCount} <(awk '{if(\$1!~/^chr/) print \'chr\'\$0; else print \$0}' '{}') $N '{}'_mapped' &&
paste {WorkDir}/*_mapped > {WorkDir}/data &&
{ProcessR} {GinkgoHomeDir}/genomes/{Genome}/original {WorkDir} {R_Args_Base} && {ReclustR} {GinkgoHomeDir}/genomes/{Genome}/original {WorkDir} {R_Args_Reclust} && {CNVCaller} {WorkDir}/SegCopy {WorkDir}/CNV1 {WorkDir}/CNV2"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    if ReclustR:
        _tgt = os.path.dirname(ReclustR) if os.path.splitext(ReclustR)[1] else ReclustR
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if GinkgoHomeDir:
        _tgt = os.path.dirname(GinkgoHomeDir) if os.path.splitext(GinkgoHomeDir)[1] else GinkgoHomeDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if WorkDir:
        _tgt = os.path.dirname(WorkDir) if os.path.splitext(WorkDir)[1] else WorkDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if NormalDir:
        _tgt = os.path.dirname(NormalDir) if os.path.splitext(NormalDir)[1] else NormalDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BinFile:
        _tgt = os.path.dirname(BinFile) if os.path.splitext(BinFile)[1] else BinFile
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if CNVCaller:
        _tgt = os.path.dirname(CNVCaller) if os.path.splitext(CNVCaller)[1] else CNVCaller
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BamToBedDir:
        _tgt = os.path.dirname(BamToBedDir) if os.path.splitext(BamToBedDir)[1] else BamToBedDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if ProcessR:
        _tgt = os.path.dirname(ProcessR) if os.path.splitext(ProcessR)[1] else ProcessR
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BinMeth:
        _tgt = os.path.dirname(BinMeth) if os.path.splitext(BinMeth)[1] else BinMeth
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BinUnsorted:
        _tgt = os.path.dirname(BinUnsorted) if os.path.splitext(BinUnsorted)[1] else BinUnsorted
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if BinCount:
        _tgt = os.path.dirname(BinCount) if os.path.splitext(BinCount)[1] else BinCount
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    if ResultBaseDir:
        _tgt = os.path.dirname(ResultBaseDir) if os.path.splitext(ResultBaseDir)[1] else ResultBaseDir
        if _tgt: os.makedirs(_tgt, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()