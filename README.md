# Genetic QC Nextflow Pipeline

A reproducible Nextflow DSL2 pipeline for genetic quality control.

The project contains two independent workflows:

- `snp_array_qc/`: QC for PLINK binary files from SNP-array or GWAS-style data.
- `wgs_wes_qc/`: QC for whole-genome or whole-exome sequencing data.

The pipeline separates variant-level QC from sample-level QC. This is important because variant-level QC can often be run per chromosome, while sample-level QC usually requires genome-wide context.

All thresholds are defined in config files and can be overridden at runtime. Individual modules and whole QC phases can be enabled or disabled with parameters.

## Workflow Design

Both workflows follow the same high-level structure:

```text
01_input_check
02_variant_level_qc
03_sample_level_qc
04_final_report
```

Variant-level QC includes filters and metrics applied to variants, sites, or genotypes. These steps can often run per chromosome:

- duplicate variant checks
- variant missingness
- HWE filtering
- MAF filtering
- site-level filtering
- genotype-level filtering
- Ti/Tv and variant count summaries
- chromosome merge or aggregation

Sample-level QC includes filters and metrics applied to individuals. These steps should be treated as final only when all autosomes are present:

- sample missingness
- heterozygosity
- sex check
- relatedness
- ancestry PCA
- contamination
- sample-level variant counts
- coverage and alignment metrics

The `sample_qc_scope` parameter controls this behavior:

| Value | Behavior |
|-------|----------|
| `auto` | Use `genome_wide` when all autosomes are present; otherwise mark sample QC as provisional. |
| `genome_wide` | Treat sample-level QC as final. |
| `provisional` | Run sample-level QC but mark it as not suitable for final filtering. |
| `skip` | Skip sample-level QC. |

## Setup

Install one execution environment:

```bash
bash setup.sh
bash setup.sh docker
bash setup.sh singularity
```

Then test the tools:

```bash
bash test_env.sh
bash test_env.sh docker
bash test_env.sh singularity
```

Each tool is reported as `PASS`, `WARN`, or `FAIL`. The full log is written to `test_results.log`.

Example output:

```text
[PASS]  PLINK 1.9              PLINK v1.90b6.21
[PASS]  PLINK2                 PLINK v2.00a5.12
[PASS]  samtools               samtools 1.18
[PASS]  bcftools               bcftools 1.18
[PASS]  picard                 3.1.1
[PASS]  FastQC                 FastQC v0.12.1
[PASS]  mosdepth               mosdepth 0.3.6
[PASS]  GATK                   4.5.0.0
[WARN]  VerifyBamID2           not found
[PASS]  Python 3               Python 3.10.12
[PASS]  R                      R version 4.3.1
Results: 13 PASS  1 WARN  0 FAIL
```

If any required tool reports `FAIL`, re-run `setup.sh`, inspect the error, or check [containers/environment.yml](containers/environment.yml).

For HPC use, after `bash setup.sh singularity`, copy `containers/genetic-qc.sif` to shared cluster storage and update `conf/singularity.config`.

## SNP Array QC

Run the SNP-array workflow:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  -profile docker
```

Run with a reference panel for ancestry PCA:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --reference_panel data/1000G/1000G_hg38 \
  --outdir results/snp_array_qc \
  -profile singularity
```

Run a chromosome test and mark sample QC as provisional:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --chroms 22 \
  --sample_qc_scope provisional \
  --outdir results/snp_array_qc_chr22
```

Skip sample-level QC:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --run_sample_qc false \
  --outdir results/snp_array_variant_only
```

Override thresholds:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --maf 0.05 \
  --hwe_p 1e-4 \
  --sample_missingness 0.05 \
  --outdir results/snp_array_qc_loose
```

## WGS / WES QC

Run WES QC from BAM input:

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

Run VCF-based variant QC by chromosome, then merge and run genome-wide sample QC:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --mode wgs \
  --chroms 1-22 \
  --outdir results/vcf_qc \
  -profile docker
```

Run a chromosome-only test:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --mode wgs \
  --chroms 22 \
  --sample_qc_scope provisional \
  --outdir results/vcf_chr22_test
```

## Params Files

Run with a params file:

```bash
nextflow run snp_array_qc/main.nf -params-file assets/example_params.yml
```

See [assets/example_params.yml](assets/example_params.yml) for an annotated example.

## Execution Profiles

Profiles can be combined, for example `-profile slurm,singularity`.

| Profile | When to use |
|---------|-------------|
| `standard` | Tools are already available on `PATH`. |
| `docker` | Local machine or cloud execution with Docker. |
| `singularity` | HPC clusters using Singularity or Apptainer. |
| `conda` | Systems without a container daemon. |
| `slurm` | SLURM scheduler; usually combine with a container profile. |
| `lsf` | IBM LSF scheduler. |
| `awsbatch` | AWS Batch. |

Edit queue and account settings in `nextflow.config`:

```groovy
slurm {
    process.queue          = 'your_queue'
    process.clusterOptions = '--account=your_account'
}
```

## Key Parameters

### Shared

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_variant_qc` | true | Enable variant-level QC phase. |
| `--run_sample_qc` | true | Enable sample-level QC phase. |
| `--chroms` | 1-22 | Chromosomes to process. Accepts `1-22`, `22`, `1,2,22`, or `all`. |
| `--sample_qc_scope` | auto | One of `auto`, `genome_wide`, `provisional`, or `skip`. |
| `--keep_intermediate` | false | Retain intermediate files. |
| `--outdir` | results | Output directory. |

### SNP Array

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bfile` | required | PLINK binary prefix without `.bed`, `.bim`, `.fam`. |
| `--pheno` | null | Optional phenotype or covariate file. |
| `--reference_panel` | null | Optional reference panel prefix for ancestry PCA. |
| `--sample_missingness` | 0.02 | Maximum missing genotype fraction per sample. |
| `--variant_missingness` | 0.02 | Maximum missing genotype fraction per variant. |
| `--hwe_p` | 1e-6 | HWE p-value cutoff. |
| `--maf` | 0.01 | Minimum minor allele frequency. |
| `--heterozygosity_sd` | 3 | Heterozygosity outlier SD threshold. |
| `--relatedness_pi_hat` | 0.1875 | Relatedness cutoff. |
| `--pca_outlier_sd` | 6 | PCA outlier SD threshold. |

### WGS / WES

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_type` | bam | One of `fastq`, `bam`, `cram`, or `vcf`. |
| `--mode` | wes | One of `wes` or `wgs`. |
| `--samplesheet` | required | CSV with `sample,file1[,file2]`. |
| `--reference_fasta` | required | Reference FASTA. A `.fai` index is expected. |
| `--target_intervals` | required for WES | Capture target BED file. |
| `--min_mean_depth_wes` | 30 | Minimum mean on-target depth. |
| `--min_mean_depth_wgs` | 20 | Minimum mean genome-wide depth. |
| `--min_target_20x_fraction` | 0.80 | Minimum fraction of targets covered at least 20x. |
| `--max_contamination` | 0.03 | Maximum allowed contamination estimate. |
| `--max_duplication_rate` | 0.20 | Maximum allowed duplicate rate. |
| `--min_gq` | 20 | Minimum genotype quality. |
| `--min_dp` | 10 | Minimum genotype depth. |
| `--variant_filter_method` | hard_filter | One of `hard_filter` or `vqsr`. |
| `--pca_maf` | 0.05 | Minimum MAF for PCA variant selection. |
| `--pca_autosomes_only` | true | Restrict PCA variant selection to autosomes. |

## Output Structure

SNP-array output:

```text
results/snp_array_qc/
|-- logs/
|-- cleaned_data/
|-- exclusion_lists/
|-- qc_tables/
|-- qc_plots/
`-- final_report.html
```

WGS/WES output:

```text
results/wgs_wes_qc/
|-- logs/
|-- fastqc/
|-- alignment_metrics/
|-- duplicate_metrics/
|-- coverage_qc/
|-- contamination/
|-- sex_check/
|-- variant_calling_qc/
|-- ancestry_pca/
|-- cleaned_data/
|-- qc_tables/
|-- qc_plots/
`-- wgs_wes_final_report.html
```

## Documentation

- [SNP Array QC Manual](docs/snp_array_qc_manual.md)
- [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md)
- [References](docs/references.md)

## References

- Anderson et al. 2010, *Nature Protocols*: GWAS QC thresholds.
- GATK Best Practices: WGS/WES variant filtering and VQSR.
- gnomAD v4.0: contamination cutoffs and sample-level QC framework.
- PLINK 1.9 and PLINK 2 documentation.
