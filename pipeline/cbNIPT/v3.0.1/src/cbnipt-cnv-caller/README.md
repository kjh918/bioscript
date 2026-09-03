# cbnipt-cnv-caller

cbNIPT manual-method CNV & variant-profile caller.

## CLI

```bash
cbnipt-cnv-caller make-bins --ReferenceFasta hg38.fa --MappabilityBW hg38.100mer.bw --OutBinFile out/Binning
cbnipt-cnv-caller call-cnv  --SeqID S1 --BamPath S1.bam --VcfFile sites.vcf.gz \
                            --AnnotatedBins out/Binning/hg38.fa.bin_100.0K.bed.gz --OutDir out/S1
```

## Pipeline (`CnvPipeline.run`)

1. `binning` – load pre-computed annotated bins
2. `calling.processing` – parallel BAM/VCF evidence extraction, LOWESS GC correction
3. `calling.normalization` – QC bin filter + Leave-One-Out sex-aware normalization (baseline CN = 2.0)
4. `calling.sex_calling` – fetal sex + sex-chromosome aneuploidy hard calls (single rule table in `config.yaml`)
5. `calling.chrom_calling` – chromosome-level ABNORMAL / SUSPICIOUS / NORMAL calls
6. `calling.segmentation` – macro CNV (PELT changepoint) + micro CNV (rolling median)
7. `clinical.clinical_diagnosis` – cross-check against the NIPT marker DB (`resources/gcx_nipt_target.tsv`)
8. `viz` – QC / overview / heatmap / clinical-threshold plots

All thresholds live in `src/cbnipt_cnv_caller/config.yaml`. Override the marker DB with
`CBNIPT_MARKER_DB=/path/to/target.tsv`.

## Tests

```bash
uv run pytest packages/cbnipt-cnv-caller
```
