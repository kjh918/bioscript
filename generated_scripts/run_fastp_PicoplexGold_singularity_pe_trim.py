#!/usr/bin/env python3
import argparse, json, re, subprocess, os, shlex
from pathlib import Path

TOKEN_PAT = re.compile(r"\[([A-Za-z0-9_]*)\]")
CMD_LINE = 'singularity exec -B [bind] [sif] fastp --thread [Threads] --in1 [RawFastqDir]/[SeqID]_R1.fastq.gz --in2 [RawFastqDir]/[SeqID]_R2.fastq.gz --out1 [out_read1] --out2 [out_read2] --json [json] --html [html] --trim_poly_g --detect_adapter_for_pe --adapter_sequence [adapter_sequence] --adapter_sequence_r2 [adapter_sequence_r2] --length_required [length_required] --average_qual [average_qual] --qualified_quality_phred [qualified_quality_phred] --trim_front1 [trim_front1] --trim_front2 [trim_front2]'
REQUIRED_KEYS = ['SeqID', 'RawFastqDir', 'TrimFastqDir', 'qcResDir']
DEFAULTS = {'Threads': '8', 'sif': '/storage/images/fastp-0.23.4.sif', 'bind': '/storage,/data', 'length_required': '100', 'average_qual': '10', 'qualified_quality_phred': '15', 'trim_front1': '14', 'trim_front2': '14', 'adapter_sequence': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'adapter_sequence_r2': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'out_read1': '[TrimFastqDir]/[SeqID].trimmed_R1.fastq.gz', 'out_read2': '[TrimFastqDir]/[SeqID].trimmed_R2.fastq.gz', 'json': '[qcResDir]/[SeqID].fastp.json', 'html': '[qcResDir]/[SeqID].fastp.html'}
OUTPUT_KEYS = ['out_read1', 'out_read2', 'json', 'html']

def render(s, ctx):
    def repl(m):
        key = m.group(1)
        return str(ctx.get(key, m.group(0))).strip()
    res = s
    for _ in range(3):
        if "[" not in res: break
        res = TOKEN_PAT.sub(repl, res)
    return " ".join(res.split())

def main():
    parser = argparse.ArgumentParser()
    for k in ['RawFastqDir', 'SeqID', 'Threads', 'TrimFastqDir', 'adapter_sequence', 'adapter_sequence_r2', 'average_qual', 'bind', 'html', 'json', 'length_required', 'out_read1', 'out_read2', 'qcResDir', 'qualified_quality_phred', 'sif', 'trim_front1', 'trim_front2']:
        parser.add_argument(f"--{k}", default=DEFAULTS.get(k, ""))
    parser.add_argument("--cwd", default=".")
    parser.add_argument("--emit-outputs", default="")

    args = parser.parse_args()
    ctx = vars(args)

    for rk in REQUIRED_KEYS:
        if not ctx.get(rk):
            print(f"Error: Missing --{rk}")
            exit(1)

    # Output 경로 렌더링
    for k in OUTPUT_KEYS:
        ctx[k] = render(ctx[k], ctx)

    final_cmd = render(CMD_LINE, ctx)
    print(f"\n[RUNNING CMD]\n{final_cmd}\n")

    os.makedirs(args.cwd, exist_ok=True)
    proc = subprocess.run(final_cmd, shell=True, cwd=args.cwd)

    if args.emit_outputs:
        out_data = {k: ctx[k] for k in OUTPUT_KEYS}
        Path(args.emit_outputs).write_text(json.dumps(out_data, indent=2))

    exit(proc.returncode)

if __name__ == "__main__":
    main()
