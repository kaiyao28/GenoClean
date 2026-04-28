# Genetic QC Pipeline

Everything you need for production-grade genetic QC, batteries included. Clone the repo, run one command, get a clean dataset and an HTML report — no manual tool installation, no custom scripts.

The pipeline covers the full QC stack for two data types:

| Workflow | Input | Use case |
|----------|-------|----------|
| `snp_array_qc/` | PLINK binary (`.bed/.bim/.fam`) | SNP arrays, GWAS datasets |
| `wgs_wes_qc/` | FASTQ, BAM/CRAM, or VCF | Whole-genome or whole-exome sequencing |

Every QC step — variant missingness, HWE, MAF, sex check, heterozygosity, relatedness, ancestry PCA, contamination, coverage, duplicate rate — runs in the right order, with the right tools, all pre-installed in a single Docker image. Each step can be toggled on or off independently, and every threshold has a sensible default that can be overridden from the command line. The pipeline resumes from where it left off if anything fails.

For SNP array data, the recommended workflow has two stages: **inspect first, then filter**. The inspection step (`inspect.nf`) computes all metric distributions on unfiltered data and produces an annotated `params_template.yaml` so you can choose appropriate thresholds for your specific dataset before anything is removed. The QC step (`main.nf`) then applies those thresholds. Both steps are fast; running them sequentially takes little more time than running QC directly, and it eliminates the risk of applying the wrong threshold to your data.

At the end you get a self-contained HTML report with the full attrition table, all metrics, and the exact thresholds used — ready to paste into a methods section.

---

## Quick Start

```bash
git clone https://github.com/kaiyao28/GeneticQC.git
cd GeneticQC
bash test_data/run_smoke_tests.sh              # both pipelines
bash test_data/run_smoke_tests.sh --test snp_array   # SNP array only
bash test_data/run_smoke_tests.sh --test wgs_wes     # WGS/WES only
```

The smoke test script checks Docker and Nextflow, pulls the image, and runs the selected workflow(s) on synthetic toy data in `test_data/`. Both tests should complete in a few minutes and write HTML reports to `results/`.

For HPC without Docker or Singularity, use `--profile manual_paths` instead (see [Setup Guide](docs/setup.md)).

For platform-specific setup (Windows/WSL, Linux, macOS, HPC), see [Setup Guide](docs/setup.md).

---

## Example: SNP Array QC

### Step 1 — Inspect your data first

```bash
nextflow run snp_array_qc/inspect.nf \
  --bfile data/raw/genotypes \
  --outdir results/inspect \
  -profile docker
```

Open `results/inspect/inspect_report.html` in a browser. It shows the full distribution of every QC metric — missingness, MAF, HWE p-values, heterozygosity, pairwise IBD, and PCA — with the default thresholds marked. The file `results/inspect/params_template.yaml` is pre-filled with all parameters and annotated with the observed statistics from your dataset:

```yaml
sample_missingness: 0.02  # observed 95th pct=0.003, max=0.018 — very clean; could tighten to 0.01
hwe_p: 1.0e-6             # observed min p (autosomes)=3.2e-12, variants below 1e-6: 125
relatedness_pi_hat: 0.1875
# IBD pairs > 0.1875 (default):  12
# IBD duplicates/MZ twins:        1
```

Edit any threshold that looks wrong for your data, then fill in optional inputs:

```yaml
reference_panel: data/1000G/1000G_hg38   # for ancestry-labelled PCA
ld_regions: data/high_ld_hg19.txt        # recommended: exclude MHC and inversions
```

### Step 2 — Run QC

```bash
nextflow run snp_array_qc/main.nf \
  -params-file results/inspect/params_template.yaml \
  -profile docker
```

If you are confident the defaults are appropriate and want to skip the inspection step:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  -profile docker
```

Full parameter reference: [SNP Array QC Manual](docs/snp_array_qc_manual.md)

---

## Example: WGS / WES QC

BAM/CRAM input (WES):

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type bam \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --target_intervals /data/reference/exome_targets.bed \
  --mode wes \
  --outdir results/wgs_wes_qc \
  -profile docker
```

VCF input (variant QC only):

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --mode wgs \
  --run_variant_qc true \
  --run_sample_qc false \
  --outdir results/vcf_qc \
  -profile docker
```

Full parameter reference: [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md)

---

## Reports

Both workflows produce a self-contained HTML report at the end of the run.

**SNP Array report** (`final_report.html`):

| Metric | Value |
|--------|-------|
| Final samples | 950 |
| Final variants | 487 203 |
| Excluded samples | 12 |
| Excluded variants | 4 891 |

Followed by a per-step attrition table showing variants and samples removed at each QC filter (duplicate check, missingness, HWE, MAF, sex check, heterozygosity, relatedness, PCA), and a thresholds table recording all parameter values used.

**WGS/WES report** (`wgs_wes_final_report.html`):

| Phase | Step | Metric | Value |
|-------|------|--------|-------|
| input_check | input_check | status | PASS |
| variant_level_qc | variant_calling_qc | ts_tv_ratio | 2.07 |
| variant_level_qc | variant_calling_qc | n_snps | 4 812 301 |
| variant_level_qc | variant_calling_qc | n_indels | 891 204 |
| variant_level_qc | merge_chromosomes | n_variants | 5 703 505 |
| sample_level_qc | coverage_qc | mean_depth | 34.2 |
| sample_level_qc | contamination | freemix | 0.009 |

Followed by thresholds and run settings. The report records which phases ran, which were skipped, and whether sample-level QC is final or provisional (provisional when fewer than all 22 autosomes were analysed).

---

## Execution Profiles

| Profile | When to use |
|---------|-------------|
| `docker` | Laptop or workstation with Docker Desktop |
| `slurm,singularity` | HPC cluster with SLURM + Apptainer/Singularity |
| `lsf,singularity` | HPC cluster with LSF + Apptainer/Singularity |
| `slurm,manual_paths` | HPC with no container engine; tools installed manually |

On HPC, always pair the scheduler profile (`slurm`, `lsf`) with the container profile (`singularity`). `-profile singularity` alone runs on the login node. Use absolute paths for `--bfile` and `--outdir` on clusters — relative paths can fail silently on compute nodes.

If no container engine is available, run `bash setup_hpc_manual.sh` first to download all tools. See [Setup Guide](docs/setup.md) for full cluster instructions.

---

## Key Thresholds

All thresholds have defaults and can be overridden on the command line:

```bash
# SNP array
--maf 0.05  --hwe_p 1e-4  --sample_missingness 0.05

# WGS/WES
--min_mean_depth_wgs 30  --max_contamination 0.02  --min_gq 30
```

---

## Documentation

- [Setup Guide](docs/setup.md)
- [SNP Array QC Manual](docs/snp_array_qc_manual.md)
- [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md)
- [References](docs/references.md)
- [Test Data](test_data/README.md)
