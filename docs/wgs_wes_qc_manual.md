# WGS / WES QC Manual

## Overview

The WGS/WES QC pipeline accepts FASTQ, BAM/CRAM, or VCF input and applies a series of sample-level
and variant-level quality control steps following GATK Best Practices and gnomAD QC methods.

WES and WGS modes share the same pipeline; the key difference is:
- **WES**: requires `--target_intervals` BED file; coverage and depth thresholds apply to on-target regions
- **WGS**: no intervals required; genome-wide coverage is assessed

---

## Quick start

### BAM input, WES mode

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type bam \
  --samplesheet assets/example_samplesheet.csv \
  --reference_fasta data/reference/GRCh38.fa \
  --target_intervals data/reference/exome_targets.bed \
  --mode wes \
  --outdir results/wgs_wes_qc \
  -profile singularity
```

### FASTQ input, WGS mode

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type fastq \
  --samplesheet samplesheet.csv \
  --reference_fasta reference.fa \
  --mode wgs \
  --outdir results/wgs_wes_qc
```

### VCF input (post-variant-calling QC only)

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta reference.fa \
  --mode wes \
  --run_fastqc false \
  --run_alignment_metrics false \
  --run_duplicate_metrics false \
  --run_coverage_qc false \
  --run_contamination false
```

---

## Samplesheet format

CSV with columns:

```
sample,file1,file2
SAMPLE001,/data/bams/SAMPLE001.bam,/data/bams/SAMPLE001.bam.bai
SAMPLE002,/data/fastq/SAMPLE002_R1.fastq.gz,/data/fastq/SAMPLE002_R2.fastq.gz
```

For single-end FASTQ, leave `file2` empty or omit.

---

## Module reference

### 1. Input check

**Purpose:** Verify files exist and are indexed. Confirm reference FASTA and target intervals (WES) are accessible.

**Fails if:** Any input file is missing or target intervals absent in WES mode.

**Output:** `input_summary.tsv`, `input_check_summary.txt`

---

### 2. FastQC

**Purpose:** Assess raw read quality, adapter content, GC content, and duplication rate.

**Tool:** FastQC

**How to disable:** `--run_fastqc false` (automatically skipped for BAM/VCF input)

**Output:** `*.fastqc.html`, `fastqc_summary.txt`

**Interpretation:**
- Per-base quality below Q30 → consider trimming
- Adapter content > 10% → check adapter configuration
- Unusual GC content → library preparation issue or contamination

---

### 3. Alignment metrics

**Purpose:** Compute mapping rate, properly paired reads, chimeric reads, and insert size distribution.

**Tools:** `samtools flagstat`, Picard `CollectAlignmentSummaryMetrics`, `CollectInsertSizeMetrics`

**How to disable:** `--run_alignment_metrics false`

**Output:** `flagstat.txt`, `alignment_summary_metrics.txt`, `insert_size_metrics.txt`, `alignment_summary.txt`

**Interpretation:**
- Mapping rate < 95% → wrong reference genome or severe contamination
- Narrow insert size (< 100 bp mean) → over-fragmented library
- High chimeric rate → library preparation artefacts

---

### 4. Duplicate metrics

**Purpose:** Estimate PCR/optical duplication rate and mark duplicates for downstream tools.

**Tool:** Picard `MarkDuplicates`

**Default threshold:** `params.max_duplication_rate = 0.20` (20%)

**How to change:**
```bash
nextflow run main.nf --max_duplication_rate 0.30
```

**How to disable:** `--run_duplicate_metrics false`

**Output:** `marked_duplicates.bam`, `duplicate_metrics.txt`, `duplicate_summary.txt`

**Interpretation:**
- WGS > 20%: unusually high, investigate input DNA amount
- WES > 50%: may be expected for very tight capture kits; adjust threshold accordingly
- Estimated library size: very low values suggest poor library complexity

---

### 5. Coverage QC

**Purpose:** Verify mean sequencing depth and fraction of bases/targets covered at ≥10x, ≥20x, ≥30x.

**Tool:** mosdepth

**Default thresholds:**
- `params.min_mean_depth_wgs = 20` (WGS)
- `params.min_mean_depth_wes = 30` (WES)
- `params.min_target_20x_fraction = 0.80` (80% of targets ≥20x)

**How to change:**
```bash
nextflow run main.nf --min_mean_depth_wes 40 --min_target_20x_fraction 0.90
```

**How to disable:** `--run_coverage_qc false`

**Output:** `coverage_summary.txt`, `coverage_metrics.txt`, `coverage_plot.png`

---

### 6. Contamination check

**Purpose:** Estimate fraction of cross-sample contamination using VerifyBamID2.

**Tool:** VerifyBamID2 (SVD method, no panel required)

**Default threshold:** `params.max_contamination = 0.03` (3%)

**How to change:**
```bash
nextflow run main.nf --max_contamination 0.02
```

**How to disable:** `--run_contamination false`

**Output:** `contamination.selfSM`, `contamination_summary.txt`

**Interpretation:**
- FREEMIX > 0.03: contaminated sample, exclude from analysis
- FREEMIX 0.01–0.03: borderline; may still be usable depending on analysis

---

### 7. Sex check

**Purpose:** Infer biological sex from X and Y chromosome coverage ratios.

**Method:**
- Compute mean depth on chr1 (autosomal reference), chrX, chrY
- Female: X/auto ≈ 1.0, Y/auto ≈ 0
- Male: X/auto ≈ 0.5, Y/auto ≈ 0.5

**How to disable:** `--run_sex_check_wgs false`

**Output:** `sex_check_results.txt`, `sex_check_summary.txt`

---

### 8. Variant calling QC

**Purpose:** Compute Ti/Tv ratio, SNP/indel counts, het/hom ratio, and per-sample variant statistics.

**Tool:** `bcftools stats`

**How to disable:** `--run_variant_calling_qc false`

**Expected Ti/Tv values:**
| Scope | Expected Ti/Tv |
|-------|---------------|
| WGS (whole genome) | 2.0 – 2.1 |
| WES (coding regions) | 2.8 – 3.0 |
| dbSNP variants only | ~2.0 |

Values outside these ranges suggest calling artefacts.

---

### 9. Variant filtering

**Purpose:** Apply site-level and genotype-level filters to remove low-quality variants.

**Method:** `params.variant_filter_method = "hard_filter"` (default) or `"vqsr"`

**Hard-filter thresholds (SNPs):**
```
QD < 2.0       QualByDepth: normalised variant quality
FS > 60.0      FisherStrand: strand bias
MQ < 40.0      RMSMappingQuality
MQRankSum < -12.5
ReadPosRankSum < -8.0
QUAL < 30
```

**Hard-filter thresholds (indels):**
```
QD < 2.0
FS > 200.0
ReadPosRankSum < -20.0
QUAL < 30
```

**Genotype filters:**
```
GQ < params.min_gq (default 20)
DP < params.min_dp (default 10)
```

**How to change any filter:**
```bash
nextflow run main.nf --snp_qd 3.0 --min_gq 30
```

**How to disable:** `--run_variant_filtering false`

**VQSR mode:** Set `--variant_filter_method vqsr` and provide training resource VCFs:
```bash
nextflow run main.nf \
  --variant_filter_method vqsr \
  --vqsr_hapmap /path/hapmap.vcf \
  --vqsr_omni /path/omni.vcf \
  --vqsr_1000g /path/1000G.vcf \
  --vqsr_dbsnp /path/dbsnp.vcf
```

---

### 10. Sample variant counts

**Purpose:** After genotype filtering, compute per-sample call rate and flag samples with low call rate.

**Default threshold:** `params.min_call_rate = 0.95` (95%)

**How to change:**
```bash
nextflow run main.nf --min_call_rate 0.90
```

**How to disable:** `--run_sample_variant_counts false`

**Output:** `sample_variant_counts.tsv`, `sample_count_outliers.txt`

---

### 11. Relatedness

**Purpose:** Detect related or duplicate sample pairs from the VCF.

**Tool:** `bcftools +relatedness2`

**Default threshold:** `params.relatedness_pi_hat = 0.1875`

**How to disable:** `--run_relatedness_wgs false`

---

### 12. Ancestry PCA

**Purpose:** PCA on high-quality, LD-pruned, common biallelic SNPs extracted from the VCF.

**Tool:** PLINK2 with `--pca`

**Default threshold:** `params.pca_outlier_sd = 6`

**How to disable:** `--run_ancestry_pca_wgs false`

---

## Output structure

```
results/wgs_wes_qc/
├── fastqc/                 # FastQC reports
├── alignment_metrics/      # samtools and Picard outputs per sample
├── duplicate_metrics/      # marked BAMs and duplication reports
├── coverage_qc/            # mosdepth coverage metrics
├── contamination/          # VerifyBamID2 results
├── sex_check/              # sex inference results
├── variant_calling_qc/     # bcftools stats outputs
├── cleaned_data/           # filtered VCFs
├── qc_tables/              # summary TSVs
├── qc_plots/               # PNG plots
├── logs/                   # pipeline logs
└── wgs_wes_final_report.html  # complete HTML report
```

---

## Common issues

| Symptom | Likely cause | Solution |
|---------|-------------|---------|
| Coverage FAIL for all samples | Wrong target intervals | Confirm BED coordinate system (0-based) matches BAM |
| Very high contamination | Sample mixture or labelling error | Investigate sample provenance |
| Ti/Tv < 1.5 | Calling artefacts or no reference | Check GATK parameters and reference genome |
| Low call rate after genotype filter | Too stringent GQ/DP | Try `--min_gq 10 --min_dp 5` for lower-coverage WES |
| VQSR fails | Cohort too small | Use `--variant_filter_method hard_filter` instead |
