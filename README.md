# nf-mutect2

Somatic SNV/indel calling with GATK Mutect2 for WGS (hg38), including optional Panel of Normals (PON), scatter/gather sharding, and filtering. Designed for SLURM but runnable locally.

## Requirements
- Nextflow >= 22.10
- Java 11+
- Conda (default) or Singularity/Docker if you adapt `nextflow.config`
- Reference: hg38 FASTA with .fai and .dict (created automatically)
- Known sites VCF for contamination (and for BQSR when enabled) (e.g., gnomAD/common)

## Sample sheet (preferred)
CSV with header:
```
subject,type,sample_id,fastq1,fastq2,cram
```
- `type`: `tumor` or `normal`
- Supply either paired FASTQs (fastq1/fastq2) or a CRAM path per row.
- Tumor-only is allowed: include a `tumor` row without a matching `normal` for the same `subject`.

Example:

```
subject,type,sample_id,fastq1,fastq2,cram

PAT1,tumor,TUMOR1,/data/TUMOR1_R1.fq.gz,/data/TUMOR1_R2.fq.gz,

PAT1,normal,NORM1,/data/NORM1_R1.fq.gz,/data/NORM1_R2.fq.gz,

PAT2,tumor,TUMOR2,,,/crams/TUMOR2.cram

```

## Quickstart (SLURM profile)
```bash
# PON build (optional, normals already CRAMed & BQSR-ready)
nextflow run mutect2_pipeline.nf -profile slurm \
  --reference /refs/hg38.fa \
  --known_sites /refs/gnomad.vcf.gz \
  --pon_crams "/normals/*.cram" \
  --scatter_count 50 \
  --outdir results/pon

# Tumor/normal calling with scatter/gather + filtering
nextflow run mutect2_pipeline.nf -profile slurm \
  --reference /refs/hg38.fa \
  --known_sites /refs/gnomad.vcf.gz \
  --sample_sheet samples.csv \
  --panel_of_normals results/pon/panel_of_normals.vcf.gz \
  --run_bqsr false \
  --scatter_count 50 \
  --outdir results/case01
```

If you start from FASTQs, swap `--crams` for `--reads "data/*_{R1,R2}.fastq.gz"`.

## Dry run (stub)
```bash
nextflow run mutect2_pipeline.nf -profile local,stub -stub-run \
  --reference /refs/hg38.fa \
  --known_sites /refs/gnomad.vcf.gz \
  --sample_sheet samples.csv \
  --outdir results/stub
```
If `nextflow` is not on `PATH`, use your local binary, for example:
`/Users/adigenova/nextflow_pmd/nextflow run mutect2_pipeline.nf ...`

## Outputs (under `--outdir`)
- `variants/filtered/*_filtered.vcf.gz` – final somatic calls
- `variants/shards/` – per-interval Mutect2 shards (for debugging)
- `pon/` – PON shards, merged normals, final `panel_of_normals.vcf.gz` (if built)
- `align/` – BQSR CRAMs (only when `--run_bqsr true`)
- `qc/fastqc`, `qc/multiqc` – QC reports
- `timeline.html`, `report.html`, `trace.txt`, `dag.svg` – execution metadata

## Notes
- `--scatter_count` defaults to 50 for WGS; reduce for small panels.
- Provide `--intervals` (BED or interval_list) to restrict calling; otherwise SplitIntervals shards the whole reference.
- Duplicate marking is assumed to be done during alignment (no extra GATK MarkDuplicates pass).
- BQSR is optional and disabled by default. Enable with `--run_bqsr true`.
- Known sites are required for contamination estimation (and for BQSR when enabled).
- Resource defaults live in `nextflow.config`; tune per cluster.
- Tumor-only mode is supported when no normal exists for a subject in the sample sheet; contamination is computed tumor-only and Mutect2 runs without `-normal`.
