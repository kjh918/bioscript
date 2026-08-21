#!/bin/bash
# ==========================================
# Tool Setup : SRAToolkit (v3.0.10)
# Profile    : SRA Download, FASTQ Conversion, Rename, and Gzip Compression
# Writer     : jhkim (2026-07-29)
# ==========================================

mkdir -p [OutDir] [TmpDir]
test -x [prefetch_bin]
test -x [fasterq_dump_bin]
