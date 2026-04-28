# SNP Array QC Manual

## Overview

The SNP array QC pipeline processes PLINK binary files (`.bed`/`.bim`/`.fam`) through a series
of established quality control steps. Each step can be enabled/disabled independently, and every
threshold is configurable without editing process code.

The module order follows Anderson et al. 2010 (Nature Protocols), which is the most-cited
reference for GWAS quality control procedures.

---

## Before you start

### Checklist

- [ ] **Nextflow installed** — `nextflow -version` should show ≥ 22.10. See [Setup Guide](setup.md).
- [ ] **Container available** — `singularity --version` on HPC, or `docker version` on a laptop.
- [ ] **Data in PLINK binary format** — three files with the same prefix: `.bed`, `.bim`, `.fam`.
- [ ] **Chromosome coding is consistent** — chromosomes must be numbered `1 2 ... 22 23`, not `chr1 chr2 ... chrX`. Check with `awk '{print $1}' genotypes.bim | sort -u`.
- [ ] **FAM file column 6 is correct** — this controls case-control mode (see below).

### Converting to PLINK binary format

If your data is in VCF format:
```bash
plink2 --vcf genotypes.vcf.gz --make-bed --out genotypes
```

If you have PLINK text format (`.ped`/`.map`):
```bash
plink --file genotypes --make-bed --out genotypes
```

### FAM file column 6 (phenotype coding)

The FAM file has six columns: `FID IID PID MID SEX PHENO`. Column 6 (PHENO) tells the pipeline whether this is a case-control study:

```
FAM001 S001 0 0 1 1   ← control (1)
FAM001 S002 0 0 2 2   ← case (2)
FAM001 S003 0 0 1 0   ← missing phenotype (0 or -9)
```

- `1` = control, `2` = case → **case-control mode**: HWE tested in controls only; differential missingness test applied
- `0` or `-9` = missing → **cohort mode**: HWE tested in all samples; differential missingness test not applied

If your FAM file uses `0`/`1` coding instead of `1`/`2`, the pipeline will not detect case-control status. Recode before running:
```bash
awk 'BEGIN{OFS="\t"} {if($6==1) $6=2; else if($6==0) $6=1; print}' genotypes.fam > genotypes_recoded.fam
```

---

## Running on a cluster

### Choose your profile

| Where | Profile |
|-------|---------|
| Laptop with Docker | `-profile docker` |
| HPC with Apptainer/Singularity + SLURM | `-profile slurm,singularity` |
| HPC with Apptainer/Singularity + LSF | `-profile lsf,singularity` |
| HPC with no container (tools installed manually) | `-profile slurm,manual_paths` |

On HPC, always combine the **scheduler profile** (`slurm`, `lsf`) with the **container profile** (`singularity`). Using `-profile singularity` alone runs everything on the login node, which is slow and may get your session terminated.

### Use absolute paths

On clusters, compute nodes may have different working directories. Use absolute paths for `--bfile` and `--outdir`:

```bash
nextflow run snp_array_qc/inspect.nf \
  --bfile /shared/data/project/genotypes \
  --outdir /shared/results/inspect \
  -profile slurm,singularity
```

Relative paths (`data/genotypes`) work on a laptop but can fail silently on cluster compute nodes.

### Submit Nextflow itself as a job

On many clusters, running long commands on the login node is not allowed. Wrap the `nextflow run` command in a job submission:

```bash
# SLURM
sbatch --mem=4G --time=24:00:00 --wrap="
  nextflow run snp_array_qc/main.nf \
    -params-file /shared/results/inspect/params_template.yaml \
    -profile slurm,singularity \
    -resume
"
```

The Nextflow head process is lightweight — it only coordinates. The actual PLINK and R work runs on separate compute nodes submitted automatically.

### Resume after a failure

Add `-resume` to any command. Nextflow caches completed steps and skips them on re-run:

```bash
nextflow run snp_array_qc/main.nf \
  -params-file results/inspect/params_template.yaml \
  -profile slurm,singularity \
  -resume
```

### Work directory

Nextflow writes intermediate files to `work/` in the directory where you run the command. On HPC, redirect this to scratch storage to avoid filling your home quota:

```bash
nextflow run ... -work-dir /scratch/$USER/nf-work
```

---

## Recommended workflow

**Always run `inspect.nf` before `main.nf`.** QC thresholds are not universal — a missingness
cutoff that is appropriate for a well-genotyped clinical array can be far too lenient or too
strict for a different platform, batch, or population. Running `inspect.nf` first takes only a
few minutes and gives you the full picture of your data before any sample or variant is removed.

```
inspect.nf  →  review  →  edit params_template.yaml  →  main.nf
```

Running `main.nf` directly with defaults is supported and works for most datasets, but it is
best reserved for cases where you already know your data quality profile.

### Step 1 — Inspect (no filtering)

```bash
nextflow run snp_array_qc/inspect.nf \
  --bfile data/raw/genotypes \
  --outdir results/inspect \
  -profile singularity
```

This takes roughly the same time as a single PLINK run. Three processes execute in parallel:

| Process | What it computes |
|---------|-----------------|
| Stats | Per-sample missingness, heterozygosity, sex check F-statistic; per-variant missingness, MAF, HWE p-values |
| Relatedness | Pairwise IBD scan; pair counts at 0.125 / 0.1875 / 0.25 / 0.5 / 0.9 thresholds |
| PCA | 10-PC PCA; variance explained per PC |

**Outputs in `results/inspect/`:**

| File | Description |
|------|-------------|
| `inspect_report.html` | Self-contained HTML — all distribution plots with default thresholds marked, plus observed stats table. Open this first. |
| `params_template.yaml` | Every parameter pre-filled at its default, annotated with observed statistics from your specific dataset |
| `qc_plots/*.png` | Individual plot files for each metric |

### Step 2 — Review

Open `inspect_report.html` in a browser. For each metric, the plot shows where the default
threshold falls relative to the actual distribution. Ask:

- **Sample / variant missingness:** Does the default (0.02) sit above the main body of the
  distribution? If most samples are below 0.005, tighten to 0.01.
- **MAF:** Is there a large spike of very rare variants you want to exclude? Or are you doing
  rare-variant analysis where you want to keep them?
- **HWE:** How many variants fall below `1e-6`? A handful is normal (genuine association signals
  or genotyping error). Hundreds may indicate batch problems or population stratification.
- **Heterozygosity:** Is the distribution roughly bell-shaped? Bimodal distribution suggests
  mixed ancestry — run PCA first and apply het checks within ancestry groups.
- **IBD:** Are there more related pairs than expected? If duplicates or MZ twins exist, decide
  whether to remove one or exclude the pair entirely.
- **PCA:** Does PC1 vs PC2 show one cloud or several? Multiple clusters indicate population
  structure — add a `--reference_panel` in main.nf to label ancestry groups.

### Step 3 — Edit the template (if needed)

The `params_template.yaml` is annotated with your observed data:

```yaml
sample_missingness: 0.02  # observed 95th pct=0.003, max=0.018 — very clean; could tighten to 0.01
variant_missingness: 0.02  # observed 95th pct=0.015, max=0.089 — default 0.02 looks appropriate
hwe_p: 1.0e-6              # observed min p (autosomes)=3.2e-12, variants below 1e-6: 125
relatedness_pi_hat: 0.1875
# IBD pairs > 0.125 (3rd degree): 18
# IBD pairs > 0.1875 (default):   12
# IBD duplicates/MZ twins:        1
```

Edit any value directly. Then fill in optional inputs if you have them:

```yaml
reference_panel: data/1000G_hg19   # for ancestry-labelled PCA plot
hapmap_info: data/hapmap_info.txt   # for HapMap population colours
ld_regions: data/high_ld_hg19.txt  # recommended — excludes MHC and chr8/17 inversions
```

### Step 4 — Run QC

```bash
nextflow run snp_array_qc/main.nf \
  -params-file results/inspect/params_template.yaml \
  -profile singularity
```

If you decide the defaults are fine as-is, you can skip editing and run directly:

```bash
nextflow run snp_array_qc/main.nf \
  -params-file results/inspect/params_template.yaml
```

Or, if you prefer not to use the template at all:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  --maf 0.05 --hwe_p 1e-4
```

---

## Module reference

### 1. Input check

**Purpose:** Validate that `.bed`, `.bim`, `.fam` exist and are non-empty. Count samples, variants, and chromosomes present.

**Command:** File existence checks + `wc -l` on `.fam` and `.bim`

**Output:**
- `input_summary.txt` — sample count, variant count, chromosome distribution, phenotype counts

**How to disable:** Cannot be disabled (always runs as the first step).

---

### 2. Imputation filter *(disabled by default)*

**Purpose:** Remove poorly imputed variants based on an imputation quality score. Only relevant for post-imputation QC — skip this step for directly genotyped SNP array data.

**Supported formats** (auto-detected from file header):

| Software | File | Quality column |
|----------|------|----------------|
| Michigan / Minimac4 | `.info` or `.info.gz` | `Rsq` |
| IMPUTE2 / IMPUTE5 | `.info` or `.info.gz` | `info` (column 5) |
| BEAGLE | `.vcf` or `.vcf.gz` | `DR2` in VCF INFO field |

**Default threshold:** `params.imputation_r2 = 0.3`
- 0.3 — Michigan's "acceptable quality" cutoff for common variants
- 0.8 — stringent; recommended for rare variant analysis

**How to enable:**
```bash
nextflow run main.nf --run_imputation_filter true \
                     --info_file results/imputed/chr_all.info.gz

# Stricter threshold
nextflow run main.nf --run_imputation_filter true \
                     --info_file data/imputed.info \
                     --imputation_r2 0.8
```

**Note on SNP ID matching:** IDs in the info file must match column 2 of the `.bim` file. Michigan uses `chr:pos:ref:alt` format — ensure PLINK files use the same convention (`plink2 --set-all-var-ids @:#:\$r:\$a`). The module prints a warning if fewer than 10% of info IDs are found in the `.bim` file.

**Output:**
- Filtered PLINK files
- `imputation_fail.txt` — variant IDs below R2 threshold
- `imputation_filter_summary.txt` — format detected, counts, ID overlap %
- `imputation_r2_plot.png` — R2 score distribution histogram

---

### 3. Duplicate variant check

**Purpose:** Remove SNPs with duplicate IDs. Arrays sometimes have redundant probes for the same variant.

**Command:** `awk` to find duplicate column-2 values in `.bim`, then `plink --exclude`

**Default behaviour:** Remove all copies of a duplicate ID (conservative).

**Output:**
- `duplicate_variants.txt` — removed variant IDs
- `duplicate_summary.txt`

**How to disable:** `--run_duplicate_check false`

---

### 4. Sample missingness

**Purpose:** Remove individuals with too many missing genotypes.

**Command:**
```bash
plink --bfile input --mind 0.02 --make-bed --out output
```

**Default threshold:** `params.sample_missingness = 0.02` (2%)

**Interpretation:** A sample with >2% missing genotypes usually reflects poor DNA quality or hybridisation failure.

**How to change:**
```bash
nextflow run main.nf --sample_missingness 0.05
```

**How to disable:** `--run_sample_missingness false`

**Output:**
- Filtered PLINK files
- `sample_callrate_removed.txt` — FID IID of removed samples
- `sample_missingness.imiss` — per-sample missingness table

---

### 5. Variant missingness

**Purpose:** Remove SNPs with too many missing calls across samples. In case-control datasets, also flags SNPs whose call rate differs significantly between cases and controls (batch-specific genotyping failures that an overall missingness filter would miss).

**Method:**
1. Compute per-locus missingness (`plink --missing`)
2. If case-control phenotype codes present: run `plink --test-missing` on the pre-filter data
3. Apply overall missingness filter + case-control exclusions in a single PLINK call

**Default thresholds:**
- `params.variant_missingness = 0.02` (2% overall missingness)
- `params.cc_miss_p = 1e-5` (differential missingness p-value; only applied when cases and controls are present)

**How to change:**
```bash
nextflow run main.nf --variant_missingness 0.05
nextflow run main.nf --cc_miss_p 1e-4   # more permissive differential test
```

**How to disable:** `--run_variant_missingness false`

**Output:**
- Filtered PLINK files (overall missingness + case-control exclusions applied)
- `variant_callrate_removed.txt` — all variants removed at this step
- `variant_missingness.lmiss` — per-variant missingness table
- `cc_miss_removed.txt` — variants failing differential missingness test (case-control only)
- `cc_miss_summary.txt` — case/control counts, threshold, number flagged
- `vmiss_plot.png` — missingness distribution histogram
- `cc_miss_plot.png` — differential missingness −log10(p) histogram (case-control only)

---

### 6. Sex check

**Purpose:** Compare genetically inferred sex (from X-chromosome F statistic) with reported sex in the FAM file.

**Command:**
```bash
plink --check-sex 0.2 0.8 --out sex_check
```

**Default thresholds:**
- `params.sex_check_f_lower_female = 0.2` — F < 0.2 → female
- `params.sex_check_f_upper_male = 0.8` — F > 0.8 → male
- Samples with 0.2 ≤ F ≤ 0.8 are flagged as ambiguous

**Interpretation:**
- Discordant: reported sex does not match genetic sex — investigate for sample swap
- Ambiguous: F statistic falls in the middle — may reflect X-chromosome anomalies or poor X coverage

**How to change:**
```bash
nextflow run main.nf --sex_check_f_lower_female 0.3 --sex_check_f_upper_male 0.7
```

**How to disable:** `--run_sex_check false`

**Output:**
- `sex_check.sexcheck` — PLINK output with F statistics
- `sex_discordant.txt` — discordant/ambiguous sample list
- `sex_check_F_stat.png` — F statistic histogram

---

### 7. Heterozygosity

**Purpose:** Flag samples with unusually high or low genome-wide heterozygosity.

**Method:** LD-pruned SNPs → PLINK `--het` → mean ± N SD filter

**Default threshold:** `params.heterozygosity_sd = 3` (±3 SD from mean)

**Interpretation:**
- Elevated heterozygosity: DNA contamination or sample mixture
- Reduced heterozygosity: inbreeding, runs of homozygosity, or DNA quality issues

**How to change:**
```bash
nextflow run main.nf --heterozygosity_sd 2
```

**How to disable:** `--run_heterozygosity false`

**Output:**
- `heterozygosity.het` — per-sample table
- `heterozygosity_outliers.txt` — outlier sample list
- `heterozygosity_plot.png` — het rate distribution with outlier cutoffs
- `miss_het_scatter.png` — missingness vs heterozygosity scatter (Anderson et al. 2010 style)

---

### 8. Relatedness

**Purpose:** Identify and recommend removal of one sample from each closely related pair.

**Method:** LD-pruned common SNPs (MAF > 5%) → `plink --genome` → greedy removal. High-LD regions (MHC, chr8/17 inversions) are excluded before pruning when `--ld_regions` is provided; these regions inflate IBD estimates and should always be excluded if a region file is available.

**Command:**
```bash
plink --genome --min 0.1875 --out relatedness
```

**Default threshold:** `params.relatedness_pi_hat = 0.1875`

PI_HAT interpretation:
| Value | Relationship |
|-------|-------------|
| ~1.0  | Monozygotic twin / duplicate |
| ~0.5  | First-degree relative |
| ~0.25 | Second-degree relative |
| ~0.125 | Third-degree relative |

**Removal strategy:** Greedy algorithm — the sample with the most related pairs is removed first.

**How to change:**
```bash
nextflow run main.nf --relatedness_pi_hat 0.125
nextflow run main.nf --ld_regions data/high_ld_regions_hg19.txt
```

**How to disable:** `--run_relatedness false`

**Output:**
- `relatedness.genome` — pairwise PI_HAT estimates
- `relatedness_remove.txt` — suggested removal list
- `relatedness_pi_hat.png` — PI_HAT distribution histogram

---

### 9. Ancestry PCA

**Purpose:** Detect population structure and flag ancestry outliers.

**Method:** LD-pruned common SNPs (MAF > 5%, geno < 1%, HWE > 0.001) → PLINK `--pca` → outlier detection (±N SD from population centre on any PC). High-LD regions are excluded when `--ld_regions` is provided.

**Default threshold:** `params.pca_outlier_sd = 6`

**Reference panel:** If `params.reference_panel` is provided, study data are merged with the reference and PCA is performed on the combined dataset. The plot colours are determined by what information is available:
1. **HapMap info provided** (`--hapmap_info`): coloured by HapMap population (CEU, CHB, YRI, etc.) with study samples in blue
2. **Reference merged, no HapMap info**: study samples vs reference panel (grey)
3. **No reference panel**: included vs ancestry outliers (red)

**How to change:**
```bash
nextflow run main.nf --pca_outlier_sd 4
nextflow run main.nf --n_pcs 20 --n_pcs_covariates 10
nextflow run main.nf --reference_panel data/1000G --hapmap_info data/hapmap_sample_info.txt
nextflow run main.nf --ld_regions data/high_ld_regions_hg19.txt
```

**How to disable:** `--run_ancestry_pca false`

**Output:**
- `pca.eigenvec` / `pca.eigenval` — all `n_pcs` computed PCs, all samples
- `pca_covariates.txt` — **analysis-ready covariate file**: study samples only, first `n_pcs_covariates` PCs (default 10), tab-separated with header. Compatible with PLINK2 (`--covar`), REGENIE (`--covar`), BOLT-LMM (`--covarFile`), PRSice (`--cov`), SAIGE
- `ancestry_outliers.txt` — FID/IID of outlier samples
- `pca_summary.txt` — variance explained per PC (PC1–PC10), outlier count, covariate dimensions
- `pca_scree.png` — variance explained per PC with cumulative line and covariate cutoff marker
- `pca_plot.png` — PC1 vs PC2 scatter

**Params:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_pcs` | 20 | PCs to compute (more = slower but allows inspecting deeper structure) |
| `n_pcs_covariates` | 10 | PCs written to `pca_covariates.txt`; must be ≤ `n_pcs` |
| `pca_outlier_sd` | 6 | SD threshold for ancestry outlier flagging |

Use the scree plot to choose `n_pcs_covariates`: look for the elbow where additional PCs explain very little variance. For typical European-ancestry cohorts, 5–10 PCs capture most structure; for multi-ancestry datasets, 10–20 may be needed.

---

### 10. HWE filter

**Purpose:** Remove SNPs that deviate significantly from Hardy-Weinberg equilibrium.

**Method:** `plink --hardy` computes HWE statistics; an awk filter then applies thresholds by chromosome. This approach is used instead of `plink --hwe` because it allows separate thresholds for autosomes and the X chromosome, and gives direct control over the test column used.

**Case-control handling:** When phenotype codes (FAM column 6: 1=control, 2=case) are detected, HWE is assessed in **controls only** (`TEST=UNAFF`). HWE deviation in cases may reflect true disease association rather than genotyping error.

**Chromosome X handling:** chrX is automatically detected at runtime.
- If chrX SNPs are present: a more permissive threshold (`hwe_p_chrx`) is applied. Males are hemizygous on X, so PLINK computes X HWE in females only — the smaller effective N makes the test noisier.
- If chrX SNPs are absent (autosome-only arrays): the autosome threshold is applied to all variants.

**Default thresholds:**
- `params.hwe_p = 1e-6` (autosomes)
- `params.hwe_p_chrx = 1e-4` (chrX, when present)

**How to change:**
```bash
nextflow run main.nf --hwe_p 1e-4          # more permissive autosomes
nextflow run main.nf --hwe_p_chrx 1e-3    # more permissive chrX
```

**How to disable:** `--run_hwe false`

**Output:**
- `hwe_stats.hwe` — per-variant HWE statistics from PLINK
- `hwe_removed_variants.txt` — IDs of removed variants
- `hwe_summary.txt` — includes `chrx_mode` field indicating whether chrX was detected
- `hwe_plot.png` — −log10(p) distribution for autosomes

---

### 11. MAF filter

**Purpose:** Remove variants with very low minor allele frequency.

**Command:**
```bash
plink --maf 0.01 --make-bed --out output
```

**Default threshold:** `params.maf = 0.01` (1%)

**When to change:**
- Rare variant analysis (burden tests, SKAT): set `--maf 0.001` or `--run_maf_filter false`
- Common variant GWAS (standard): keep 0.01
- Stringent common-variant analysis: `--maf 0.05`

**How to disable:** `--run_maf_filter false`

**Output:**
- `maf_removed_variants.txt` — IDs of variants removed
- `maf_summary.txt` — counts before/after filter
- `maf_plot.png` — MAF distribution histogram with threshold line (pre-filter data; shows what is removed)

---

## Workflow guide: pre-imputation vs post-imputation QC

SNP array data goes through QC in two distinct phases. The modules above serve both; which ones you enable depends on the phase.

### Phase 1 — Pre-imputation QC (directly genotyped SNP array data)

Goal: clean the array data before submitting to an imputation server.

Recommended steps (run in order):

1. Input check (always on)
2. Duplicate variant check (`--run_duplicate_check true`)
3. Sample missingness (`--sample_missingness 0.02`)
4. Variant missingness (`--variant_missingness 0.02`)
5. Sex check (`--run_sex_check true`)
6. Heterozygosity (`--run_heterozygosity true`)
7. Relatedness (`--run_relatedness true`)
8. Ancestry PCA (`--run_ancestry_pca true`)
9. HWE filter (`--hwe_p 1e-6`)
10. MAF filter (`--maf 0.01`)

Leave `--run_imputation_filter false` (the default).

Before submitting to Michigan Imputation Server or similar, also apply: strand alignment (e.g. HRC Strand Check), chromosome and position update to GRCh37/38, and removal of non-polymorphic variants.

### Phase 2 — Post-imputation QC (imputed dosage data)

Goal: filter imputed variants by quality before GWAS or PRS.

The imputation server returns dosage files (`.vcf.gz` or `.dose.vcf.gz`) with per-variant quality scores. Convert to PLINK binary format first (`plink2 --vcf ... --make-pgen --dosage-erase-dosage` or keep dosages in PLINK2 format), then run:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile results/imputed/chr_all \
  --run_imputation_filter true \
  --info_file results/imputed/chr_all.info.gz \
  --imputation_r2 0.3 \
  --run_duplicate_check false \
  --run_sex_check false \
  --run_heterozygosity false \
  --run_relatedness false \
  --run_ancestry_pca false \
  --hwe_p 1e-6 \
  --maf 0.01
```

Rationale for disabled steps: sample-level QC (sex check, het, relatedness, ancestry) was already done in Phase 1 on the directly genotyped data. Re-running on imputed data is redundant. The imputation filter replaces variant missingness as the primary variant-quality gate.

**ID format note:** Michigan imputed VCFs use `chr:pos:ref:alt` IDs. Ensure the PLINK `.bim` file uses the same convention:
```bash
plink2 --pfile ... --set-all-var-ids @:#:\$r:\$a --make-bed --out chr_all
```

---

## Output structure

```
results/snp_array_qc/
├── cleaned_data/              # final QC-passed PLINK files
├── exclusion_lists/           # all_excluded_samples.txt, all_excluded_variants.txt
├── qc_tables/                 # per-step summary .txt files and PLINK raw outputs
│   ├── sample_missingness.imiss
│   ├── heterozygosity.het
│   ├── sex_check.sexcheck
│   ├── variant_missingness.lmiss
│   ├── relatedness.genome
│   ├── cc_miss_removed.txt    # case-control only
│   └── ...
├── qc_plots/                  # PNG plots (one per QC step)
├── final_report.html          # self-contained HTML: attrition table, plots, per-sample table
├── qc_report.pdf              # multi-page PDF: title page + one plot per page + batch summary
├── qc_attrition_table.tsv     # machine-readable step-by-step counts
├── qc_thresholds.tsv          # all parameters used in this run
├── qc_per_sample.tsv          # per-sample: F_MISS, HET_RATE, sex check, QC flags, batch
├── qc_per_batch.tsv           # per-batch: sample count, removed, pass rate
└── pca_covariates.txt         # GWAS/PRS-ready: study samples × first N PCs (tab-separated)
```

---

## After QC: next steps

When `main.nf` finishes, the terminal prints the location of every output file. Key outputs:

### Cleaned dataset

`cleaned_data/` contains the QC-passed PLINK binary files, ready for association analysis:
```
results/snp_array_qc/cleaned_data/genotypes_sample_qc_pass.bed
results/snp_array_qc/cleaned_data/genotypes_sample_qc_pass.bim
results/snp_array_qc/cleaned_data/genotypes_sample_qc_pass.fam
```

### PCA covariates

`pca_covariates.txt` contains the first N principal components for every QC-passed sample. Pass it directly to your analysis tool:

| Tool | Flag |
|------|------|
| PLINK2 | `--covar pca_covariates.txt` |
| REGENIE | `--covar pca_covariates.txt` |
| BOLT-LMM | `--covarFile pca_covariates.txt` |
| PRSice | `--cov pca_covariates.txt` |
| SAIGE | `covFile = "pca_covariates.txt"` |

The file contains only study-cohort samples (reference panel samples excluded if a reference was merged during PCA). Change the number of PCs written with `--n_pcs_covariates`.

### Review the report

Open `final_report.html` in a browser. It shows the full attrition table (samples and variants removed at each step), all plots, the per-sample QC table, and the exact thresholds used. If any step removed far more samples or variants than expected, re-run `inspect.nf` and adjust the relevant threshold.

### Exclusion lists

`exclusion_lists/all_excluded_samples.txt` and `all_excluded_variants.txt` contain every ID removed. Keep these files — they are needed if you need to apply the same QC to a replication cohort.

### Pheno file format

If you use an external phenotype file (`--pheno`), it must be tab-separated with header:

```
FID     IID     PHENO
FAM001  S001    1
FAM001  S002    2
FAM001  S003    NA
```

- Column 1: FID, Column 2: IID, Column 3+: phenotype(s)
- Case-control: `1` = control, `2` = case, `NA` or `-9` = missing
- Quantitative: any numeric value; `NA` = missing
- Additional columns are allowed and passed through to PLINK as covariates

---

## Common issues

| Symptom | Likely cause | Solution |
|---------|-------------|---------|
| Many samples removed at sex check | Missing chrX SNPs, or chrX coded as 23 | Run `inspect.nf` and check the sex check F-stat plot; if bimodal at 0 and 1 the data is fine |
| cc_miss_p filter not applied | FAM col 6 not in 1/2 coding | Check `inspect_stats.tsv` — if n_cc_cases=0 cohort mode is active; recode FAM or pass `--pheno` |
| Zero variants after HWE | Threshold too stringent | Try `--hwe_p 1e-4`; check `inspect_hwe.png` to see the distribution first |
| PCA shows two distinct clouds | Population stratification | Add `--reference_panel` and filter by ancestry; inspect PCA plot first |
| Heterozygosity outliers are bimodal | Mixed ancestry | Run PCA first; apply het filter within ancestry groups separately |
| Imputation filter removes nothing | SNP ID format mismatch | Check overlap % in `imputation_filter_summary.txt`; align IDs with `plink2 --set-all-var-ids @:#:\$r:\$a` |
| Imputation filter format not detected | Unexpected header | Supported: Michigan `.info`, IMPUTE2 `.info`, BEAGLE VCF with `DR2=` |
| Pipeline fails on compute nodes but not login node | Relative paths | Use absolute paths for `--bfile` and `--outdir` |
| `--chroms` format error | Using `chr1-22` instead of `1-22` | Use numeric coding: `1-22`, `22`, `1,2,22`, or `all` |
