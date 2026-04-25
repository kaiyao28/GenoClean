# Genetic QC Nextflow Pipeline

A reproducible Nextflow DSL2 pipeline for genetic quality control, with two
independent workflows:

- **`snp_array_qc/`** — QC for PLINK binary files from SNP arrays (GWAS-style)
- **`wgs_wes_qc/`** — QC for whole-genome or whole-exome sequencing data

All thresholds are stored in config files and can be overridden at runtime
without editing any code. Every QC module can be independently enabled or disabled.

---

## Setup

**Step 1 — install** (pick one):

```bash
bash setup.sh              # Mamba/Conda environment  (works everywhere)
bash setup.sh docker       # Docker image              (local / cloud)
bash setup.sh singularity  # Docker → Singularity SIF  (HPC clusters)
```

**Step 2 — test all tools** and review the log:

```bash
bash test_env.sh               # tools on PATH  (after conda setup)
bash test_env.sh docker        # tools inside Docker image
bash test_env.sh singularity   # tools inside Singularity SIF
```

Each tool is reported as `PASS`, `WARN` (optional — pipeline still runs), or `FAIL`
(required — must fix before running). A full log is written to `test_results.log`.

```
[PASS]  PLINK 1.9              PLINK v1.90b6.21
[PASS]  PLINK2                 PLINK v2.00a5.12
[PASS]  samtools               samtools 1.18
[PASS]  bcftools               bcftools 1.18
[PASS]  picard                 3.1.1
[PASS]  FastQC                 FastQC v0.12.1
[PASS]  mosdepth               mosdepth 0.3.6
[PASS]  GATK                   4.5.0.0
[WARN]  VerifyBamID2           not found (optional — contamination uses GATK fallback)
[PASS]  Python 3               Python 3.10.12
[PASS]  R                      R version 4.3.1
Results: 13 PASS  1 WARN  0 FAIL
```

If any `FAIL` lines appear, re-run `setup.sh`, check the error, or inspect
[containers/environment.yml](containers/environment.yml) for version constraints.

> **HPC note:** after `bash setup.sh singularity`, copy
> `containers/genetic-qc.sif` to shared cluster storage and update the path
> in `conf/singularity.config`.

---

## SNP Array QC — running the pipeline

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  -profile docker
```

With a 1000 Genomes reference panel for ancestry PCA:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --reference_panel data/1000G/1000G_hg38 \
  --outdir results/snp_array_qc \
  -profile singularity
```

Override thresholds:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --maf 0.05 --hwe_p 1e-4 --sample_missingness 0.05 \
  --outdir results/snp_array_qc_loose
```

Skip individual modules:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --run_ancestry_pca false --run_maf_filter false \
  --outdir results/snp_array_qc_no_pca
```

Resume an interrupted run (completed steps are not re-run):

```bash
nextflow run snp_array_qc/main.nf --bfile data/raw/genotypes --outdir results/snp_array_qc -resume
```

---

## WGS / WES QC — running the pipeline

BAM input, WES mode:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type bam \
  --samplesheet assets/example_samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --target_intervals /data/reference/exome_targets.bed \
  --mode wes \
  --outdir results/wgs_wes_qc \
  -profile singularity
```

VCF input — variant-level QC only:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --mode wes \
  --run_fastqc false --run_alignment_metrics false \
  --run_duplicate_metrics false --run_coverage_qc false \
  --run_contamination false \
  --outdir results/vcf_qc
```

SLURM cluster with Singularity:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile /shared/data/genotypes \
  --outdir /shared/results/snp_qc \
  -profile slurm,singularity -resume
```

---

## Using a params file

```bash
nextflow run snp_array_qc/main.nf -params-file assets/example_params.yml
```

See [assets/example_params.yml](assets/example_params.yml) for a fully annotated example.

---

## Execution profiles

Profiles can be combined: `-profile slurm,singularity`

| Profile | When to use |
|---------|-------------|
| `standard` | tools already on PATH |
| `docker` | local machine or cloud |
| `singularity` | HPC clusters (Singularity / Apptainer) |
| `conda` | no container daemon available |
| `slurm` | SLURM scheduler — combine with a container profile |
| `lsf` | IBM LSF scheduler |
| `awsbatch` | AWS Batch |

Edit queue name and account in `nextflow.config`:

```groovy
slurm {
    process.queue          = 'your_queue'
    process.clusterOptions = '--account=your_account'
}
```

---

## Parameter reference

### SNP array

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bfile` | required | PLINK binary prefix (no extension) |
| `--pheno` | null | phenotype / covariate TSV |
| `--reference_panel` | null | 1000G prefix for ancestry PCA |
| `--sample_missingness` | 0.02 | max fraction missing per sample |
| `--variant_missingness` | 0.02 | max fraction missing per variant |
| `--hwe_p` | 1e-6 | HWE p-value cutoff |
| `--maf` | 0.01 | minimum minor allele frequency |
| `--heterozygosity_sd` | 3 | het outlier SD threshold |
| `--relatedness_pi_hat` | 0.1875 | PI_HAT cutoff for related pairs |
| `--pca_outlier_sd` | 6 | SD threshold for ancestry outliers |
| `--run_sex_check` | true | enable/disable sex check |
| `--run_ancestry_pca` | true | enable/disable PCA |
| `--keep_intermediate` | false | retain per-step PLINK files |
| `--outdir` | results | output directory |

### WGS / WES

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_type` | bam | fastq / bam / cram / vcf |
| `--mode` | wes | wes / wgs |
| `--samplesheet` | required | CSV: sample,file1[,file2] |
| `--reference_fasta` | required | reference FASTA (.fai required) |
| `--target_intervals` | required (WES) | capture BED file |
| `--min_mean_depth_wes` | 30 | min mean on-target depth |
| `--min_mean_depth_wgs` | 20 | min mean genome-wide depth |
| `--min_target_20x_fraction` | 0.80 | fraction of targets ≥ 20x |
| `--max_contamination` | 0.03 | max FREEMIX contamination |
| `--max_duplication_rate` | 0.20 | max PCR duplication rate |
| `--min_gq` | 20 | min genotype quality per call |
| `--min_dp` | 10 | min read depth per genotype |
| `--variant_filter_method` | hard_filter | hard_filter or vqsr |
| `--pca_maf` | 0.05 | min MAF for PCA variant selection |
| `--pca_autosomes_only` | true | restrict PCA to autosomes |
| `--outdir` | results | output directory |

---

## Output structure

### SNP array

```
results/snp_array_qc/
├── logs/
│   ├── versions.txt            # tool versions used in this run
│   └── input_summary.txt       # initial sample and variant counts
├── cleaned_data/               # final QC-passed PLINK files (.bed/.bim/.fam)
├── exclusion_lists/            # all_excluded_samples.txt / all_excluded_variants.txt
├── qc_tables/                  # per-step summary tables
├── qc_plots/                   # PNG plots (missingness, het, PCA, relatedness)
└── final_report.html           # HTML report with attrition table and embedded plots
```

### WGS / WES

```
results/wgs_wes_qc/
├── logs/versions.txt
├── fastqc/                     # FastQC HTML reports (FASTQ input only)
├── alignment_metrics/          # samtools + Picard metrics per sample
├── duplicate_metrics/          # marked BAMs and duplication rate tables
├── coverage_qc/                # mosdepth depth summaries and plots
├── contamination/              # VerifyBamID2 results
├── sex_check/                  # inferred vs reported sex
├── variant_calling_qc/         # bcftools stats, Ti/Tv tables
├── ancestry_pca/               # PC scores, eigenvalues, outlier lists, plots
├── cleaned_data/               # filtered VCF (genotype-filtered, PASS only)
├── qc_tables/
├── qc_plots/
└── wgs_wes_final_report.html
```

---

## Documentation

- [SNP Array QC Manual](docs/snp_array_qc_manual.md) — module-by-module reference
- [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md) — module-by-module reference
- [References](docs/references.md) — Anderson 2010, GATK, gnomAD, and all cited tools

---

## References

Anderson et al. 2010, *Nature Protocols* — GWAS QC thresholds (SNP array defaults).
GATK Best Practices — WGS/WES variant calling, hard-filtering, and VQSR.
gnomAD v4.0 — contamination cutoffs and sample-level QC framework.
#   G e n e t i c Q C  
 