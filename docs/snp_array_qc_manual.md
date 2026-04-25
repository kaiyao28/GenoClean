# SNP Array QC Manual

## Overview

The SNP array QC pipeline processes PLINK binary files (`.bed`/`.bim`/`.fam`) through a series
of established quality control steps. Each step can be enabled/disabled independently, and every
threshold is configurable without editing process code.

The module order follows Anderson et al. 2010 (Nature Protocols), which is the most-cited
reference for GWAS quality control procedures.

---

## Quick start

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  -profile singularity
```

Override a threshold:

```bash
nextflow run snp_array_qc/main.nf --bfile ... --maf 0.05
```

Skip a module:

```bash
nextflow run snp_array_qc/main.nf --bfile ... --run_ancestry_pca false
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

### 2. Duplicate variant check

**Purpose:** Remove SNPs with duplicate IDs. Arrays sometimes have redundant probes for the same variant.

**Command:** `awk` to find duplicate column-2 values in `.bim`, then `plink --exclude`

**Default behaviour:** Remove all copies of a duplicate ID (conservative).

**Output:**
- `duplicate_variants.txt` — removed variant IDs
- `duplicate_summary.txt`

**How to disable:** `--run_duplicate_check false`

---

### 3. Sample missingness

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

### 4. Variant missingness

**Purpose:** Remove SNPs with too many missing calls across samples.

**Command:**
```bash
plink --bfile input --geno 0.02 --make-bed --out output
```

**Default threshold:** `params.variant_missingness = 0.02` (2%)

**How to change:**
```bash
nextflow run main.nf --variant_missingness 0.05
```

**How to disable:** `--run_variant_missingness false`

**Output:**
- Filtered PLINK files
- `variant_callrate_removed.txt`
- `variant_missingness.lmiss`

---

### 5. Sex check

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

### 6. Heterozygosity

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
- `heterozygosity_plot.png`

---

### 7. Relatedness

**Purpose:** Identify and recommend removal of one sample from each closely related pair.

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
```

**How to disable:** `--run_relatedness false`

**Output:**
- `relatedness.genome` — pairwise PI_HAT estimates
- `relatedness_remove.txt` — suggested removal list

---

### 8. Ancestry PCA

**Purpose:** Detect population structure and flag ancestry outliers.

**Method:** LD-pruned SNPs → PLINK `--pca` → outlier detection (±N SD from population centre)

**Default threshold:** `params.pca_outlier_sd = 6`

**Reference panel:** If `params.reference_panel` is provided, study data are merged with the reference and projected PCs are computed.

**How to change:**
```bash
nextflow run main.nf --pca_outlier_sd 4 --reference_panel data/1000G
```

**How to disable:** `--run_ancestry_pca false`

**Output:**
- `pca.eigenvec` / `pca.eigenval`
- `ancestry_outliers.txt`
- `pca_plot.png`

---

### 9. HWE filter

**Purpose:** Remove SNPs that deviate significantly from Hardy-Weinberg equilibrium.

**Command:**
```bash
plink --hwe 1e-6 --make-bed --out output
```

**Important:** In case-control studies, HWE is assessed in **controls only**. HWE deviation in cases can reflect true disease association rather than error. PLINK handles this automatically when phenotype codes (1/2) are present in the FAM file.

**Default threshold:** `params.hwe_p = 1e-6`

**How to change:**
```bash
nextflow run main.nf --hwe_p 1e-4   # more permissive
```

**How to disable:** `--run_hwe false`

---

### 10. MAF filter

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

---

## Output structure

```
results/snp_array_qc/
├── cleaned_data/           # final QC-passed PLINK files
├── exclusion_lists/        # all sample/variant exclusion lists
├── qc_tables/              # per-step summary tables and PLINK outputs
├── qc_plots/               # PNG plots
├── logs/                   # PLINK log files and pipeline logs
└── final_report.html       # complete QC report with attrition table
```

---

## Common issues

| Symptom | Likely cause | Solution |
|---------|-------------|---------|
| Many samples removed at sex check | Missing chrX SNPs | Check `.bim` for chrX variants |
| Zero variants after HWE | HWE threshold too stringent | Try `--hwe_p 1e-4` |
| PCA shows two distinct clouds | Population stratification | Add `--reference_panel` and filter by ancestry |
| Heterozygosity outliers are bimodal | Mixed ancestry dataset | Run ancestry PCA first, then het check within ancestry groups |
