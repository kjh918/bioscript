import os
import subprocess
from glob import glob
from concurrent.futures import ProcessPoolExecutor, as_completed


# ============================================================
# Config
# ============================================================

PATH = (
    #"/storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/WGS.cfDNA.ControlSamples/"
    #"/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA837410/"
    "/storage/home/jhkim/workspace/projects/NIPT/260729-GCX-NIPT-PieplineSetup/Results/PRJNA385180/AnalysisUse/"
)

PYTHON = "/storage/home/jhkim/Apps/Python-3.11.13/python"
RUN_PY = "run.py"

BW = (
    "/storage/home/jhkim/Projects/NIPT/"
    "GCX-NIPT-260121/Resources/reference/hg38.100mer.bw"
)

FASTA = (
    "/storage/references_and_index/hg38/fasta/"
    "cbNIPT/hg38.fa"
)

MARKER_BED = (
    "/storage/home/jhkim/workspace/projects/NIPT/"
    "260729-GCX-NIPT-PieplineSetup/"
    "Resources/local/refs/merged_marker.bed"
)

OUT_BASE = (
    "/storage/home/jhkim/workspace/projects/NIPT/"
    "260729-GCX-NIPT-PieplineSetup/"
    "Analysis/WPS.cnv"
)

# 한 샘플 내부 병렬 처리 수
JOBS_PER_SAMPLE = 8

# 동시에 실행할 샘플 수
# 예: 4 samples x 8 jobs = 최대 약 32 workersq
MAX_SAMPLES = 6


# ============================================================
# Worker
# ============================================================

def run_sample(bam_path):
    bam_id = os.path.basename(bam_path).split(".")[0]

    outdir = os.path.join(
        OUT_BASE,
        bam_id,
    )

    os.makedirs(
        outdir,
        exist_ok=True,
    )

    cmd = [
        PYTHON,
        RUN_PY,
        "all",

        "--bam",
        bam_path,

        "--bw",
        BW,

        "--fasta",
        FASTA,

        "--outdir",
        outdir,

        "--sample",
        bam_id,

        "--mode",
        "L",
        "S",

        "--marker-bed",
        MARKER_BED,

        "--jobs",
        str(JOBS_PER_SAMPLE),
    ]

    print(
        f"[START] {bam_id}\n"
        f"  BAM: {bam_path}\n"
        f"  OUT: {outdir}"
    )

    # 샘플별 로그 파일
    log_path = os.path.join(
        outdir,
        f"{bam_id}.pipeline.log",
    )

    with open(log_path, "w") as log_file:

        result = subprocess.run(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode == 0:
        return {
            "sample": bam_id,
            "status": "SUCCESS",
            "returncode": result.returncode,
            "log": log_path,
        }

    return {
        "sample": bam_id,
        "status": "FAILED",
        "returncode": result.returncode,
        "log": log_path,
    }


# ============================================================
# Main
# ============================================================

def main():

    bam_file_list = sorted(
        glob(
            f"{PATH}/*/bam/*.recal.filtered*.bam"
        )
    )

    print(
        f"Found {len(bam_file_list)} BAM files"
    )

    if not bam_file_list:
        print("No BAM files found.")
        return

    print(
        f"Parallel samples : {MAX_SAMPLES}"
    )

    print(
        f"Jobs per sample  : {JOBS_PER_SAMPLE}"
    )

    print(
        f"Max jobs         : "
        f"{MAX_SAMPLES * JOBS_PER_SAMPLE}"
    )

    print("=" * 70)

    results = []

    with ProcessPoolExecutor(
        max_workers=MAX_SAMPLES
    ) as executor:

        future_to_bam = {
            executor.submit(
                run_sample,
                bam_path,
            ): bam_path

            for bam_path
            in bam_file_list
        }

        for future in as_completed(
            future_to_bam
        ):

            bam_path = future_to_bam[
                future
            ]

            bam_id = os.path.basename(
                bam_path
            ).split(".")[0]

            try:

                result = future.result()

                results.append(
                    result
                )

                if (
                    result["status"]
                    == "SUCCESS"
                ):

                    print(
                        f"[DONE] "
                        f"{bam_id}"
                    )

                else:

                    print(
                        f"[FAILED] "
                        f"{bam_id} "
                        f"returncode="
                        f"{result['returncode']}"
                    )

                    print(
                        f"         log: "
                        f"{result['log']}"
                    )

            except Exception as e:

                print(
                    f"[ERROR] "
                    f"{bam_id}: {e}"
                )

                results.append(
                    {
                        "sample": bam_id,
                        "status": "ERROR",
                        "returncode": None,
                        "log": None,
                    }
                )

    # ========================================================
    # Summary
    # ========================================================

    print("\n")
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)

    success = [
        r
        for r in results
        if r["status"] == "SUCCESS"
    ]

    failed = [
        r
        for r in results
        if r["status"] != "SUCCESS"
    ]

    print(
        f"SUCCESS : "
        f"{len(success)}"
    )

    print(
        f"FAILED  : "
        f"{len(failed)}"
    )

    if failed:

        print("\nFailed samples:")

        for r in failed:

            print(
                f"  {r['sample']} "
                f"[{r['status']}] "
                f"log={r['log']}"
            )


if __name__ == "__main__":
    main()