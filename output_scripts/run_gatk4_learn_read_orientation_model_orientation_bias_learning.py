#!/usr/bin/env python3
# [METADATA]
# TOOL_NAME = gatk4_learn_read_orientation_model
# VERSION = 4.4.0.0
# THREADS = 1
# PROFILE = orientation_bias_learning

"""
Tool: gatk4_learn_read_orientation_model (4.4.0.0)
Profile: orientation_bias_learning
"""
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser(description="gatk4_learn_read_orientation_model Analysis Runner")
    
    # --- [Argument Parsing] ---
    parser.add_argument('--SeqID', required=True, default='', help='No description (Default: )')
    parser.add_argument('--qcResDir', required=True, default='', help='f1r2 파일이 있고 모델을 저장할 디렉토리 (Default: )')
    parser.add_argument('--InputSuffix', required=False, default='', help='Input f1r2 suffix (with dot if exists) (Default: )')
    parser.add_argument('--OutputSuffix', required=False, default='', help='Output model suffix (Default: )')
    parser.add_argument('--singularity_bin', required=False, default='singularity', help='No description (Default: singularity)')
    parser.add_argument('--gatk4_sif', required=False, default='/storage/images/gatk-4.4.0.0.sif', help='No description (Default: /storage/images/gatk-4.4.0.0.sif)')
    parser.add_argument('--bind', required=False, default='/storage,/data', help='No description (Default: /storage,/data)')
    parser.add_argument('--xmx_mb', required=False, default='8192', help='8GB memory from original script (Default: 8192)')
    parser.add_argument('--Threads', required=False, default='14', help='ParallelGCThreads (Default: 14)')
    
    args = parser.parse_args()

    # --- [Variable Declarations: Key = Value] ---
    SeqID = args.SeqID
    qcResDir = args.qcResDir
    InputSuffix = args.InputSuffix
    OutputSuffix = args.OutputSuffix
    singularity_bin = args.singularity_bin
    gatk4_sif = args.gatk4_sif
    bind = args.bind
    xmx_mb = args.xmx_mb
    Threads = args.Threads

    # --- [Output Paths] ---
    read_orientation_model = f"{qcResDir}/{SeqID}.{OutputSuffix}.read-orientation-model.tar.gz"

    # --- [Command Execution] ---
    cmd = f"{singularity_bin} exec -B {bind} {gatk4_sif} gatk LearnReadOrientationModel --java-options '-XX:ParallelGCThreads={Threads} -Xmx{xmx_mb}m' -I {qcResDir}/{SeqID}.{InputSuffix}f1r2.tar.gz -O {read_orientation_model}"
    
    print(f"\\n[RUNNING]\\n{cmd}\\n")
    
    os.makedirs(os.path.dirname(qcResDir) if '.' in os.path.basename(qcResDir) else qcResDir, exist_ok=True)
    
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    main()