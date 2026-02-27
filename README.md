# nf-mutect2

Somatic SNV/indel calling with GATK Mutect2 for WGS (hg38), including optional Panel of Normals (PON), scatter/gather sharding, and filtering. Designed for SLURM but runnable locally.

## Requirements
- Nextflow >= 22.10
- Java 11+
- Conda (default) or Singularity/Docker if you adapt `nextflow.config`
- Reference: hg38 FASTA with .fai and .dict (created automatically)
- FASTQ mode expects BWA-MEM2 reference index files next to the FASTA (`.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`)
- Known sites VCF for contamination (and for BQSR when enabled) (e.g., gnomAD/common)
- Germline resource VCF for Mutect2 (e.g., `af-only-gnomad.vcf.gz`)

## Sample sheet (preferred)
CSV with header:
```
subject,type,sample_id,fastq1,fastq2,cram
```
- `type`: `tumor` or `normal`
- Supply either paired FASTQs (fastq1/fastq2) or an aligned file path per row (`.cram` or `.bam` in the `cram` column).
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
  --germline_resource /refs/af-only-gnomad.vcf.gz \
  --pon_crams "/normals/*.cram" \
  --scatter_count 50 \
  --outdir results/pon

# Tumor/normal calling with scatter/gather + filtering
nextflow run mutect2_pipeline.nf -profile slurm \
  --reference /refs/hg38.fa \
  --known_sites /refs/gnomad.vcf.gz \
  --germline_resource /refs/af-only-gnomad.vcf.gz \
  --sample_sheet samples.csv \
  --panel_of_normals results/pon/panel_of_normals.vcf.gz \
  --run_bqsr false \
  --scatter_count 50 \
  --outdir results/case01
```

If you start from FASTQs, swap `--crams` for `--reads "data/*_{R1,R2}.fastq.gz"`.
For aligned inputs, index files must already exist (`.cram.crai` for CRAM, `.bam.bai` or `.bai` for BAM); the pipeline validates and uses them as-is.

## Dry run (stub)
```bash
# Run using stub fixtures committed in this repo
tests/stub/run_stub.sh

# Or run manually
nextflow run mutect2_pipeline.nf -profile local,stub -stub-run \
  --reference tests/stub/reference.fa \
  --known_sites tests/stub/known_sites.vcf.gz \
  --germline_resource tests/stub/known_sites.vcf.gz \
  --sample_sheet tests/stub/samples.csv \
  --outdir tests/stub/results
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
- `--germline_resource` is required for all Mutect2 calls.
- Resource defaults live in `nextflow.config`; tune per cluster.
- Tumor-only mode is supported when no normal exists for a subject in the sample sheet; contamination is computed tumor-only and Mutect2 runs without `-normal`.
